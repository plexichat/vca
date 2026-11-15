/*
 * uk_vector_router.c
 *
 * High-Performance Preprocessed 8D Vector-Corridor Pathfinding Engine
 * designed for Geofabrik Whole-UK scale road networks with SQLite3 storage,
 * Structure-of-Arrays (SoA) contiguous memory, and raw pointer arithmetic.
 *
 * Preprocessing Pipeline (O(V + E)):
 *   1. Conformal Transverse Mercator projection (central meridian -2.0 deg)
 *   2. Dead-end 1-core cul-de-sac peeling
 *   3. Degree-2 chain contraction into exact macro-vectors + tortuosity
 *   4. Connected components labeling (O(1) island reachability rejection)
 *   5. Contiguous SoA packing + Reverse CSR + 128x128 spatial snap grid
 *   6. Fast binary (.vgraph) cache serialization
 */

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include "sqlite3.h"
#include "zlib.h"
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#if defined(__SSE__)
#include <xmmintrin.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline float vmr_inv_sqrt(float x) {
    if (x <= 0.0f) {
        return 0.0f;
    }
#if defined(__SSE__)
    float y = _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(x)));
    y *= 0.5f * (3.0f - x * y * y);
    return y;
#else
    return 1.0f / sqrtf(x);
#endif
}

static inline float vmr_sqrt_fast(float x) {
    if (x <= 0.0f) {
        return 0.0f;
    }
    float y = vmr_inv_sqrt(x);
    return (y > 0.0f) ? (x * y) : sqrtf(x);
}

/* ========================================================================= */
/* SECTION 1: CONSTANTS, TYPES & LOGGING                                     */
/* ========================================================================= */

#define UVR_VERSION      "2.0.0"
#define VGRAPH_MAGIC     0x564D5235U /* "VMR5" */
#define VGRAPH_FORMAT_V  8U

#define VMR_INVALID      UINT32_MAX
#define EARTH_RADIUS_M   6371008.8
#define DEG_TO_RAD       (M_PI / 180.0)
#define RAD_TO_DEG       (180.0 / M_PI)

/* OSGB-centered Conformal Transverse Mercator constants for Great Britain */
#define UK_TM_LON0_DEG   (-2.0)
#define UK_TM_LAT0_DEG   (54.0)
#define UK_TM_K0         0.9996012717

#define SPATIAL_GRID_DIM 128
#define EPS_F            1e-6f

typedef enum {
    LOG_ERR_L  = 0,
    LOG_WARN_L = 1,
    LOG_INFO_L = 2,
    LOG_DBG_L  = 3
} LogLvl;

static LogLvl g_log_level = LOG_INFO_L;

static void uvr_log(LogLvl lvl, const char *tag, const char *fmt, ...) {
    if (lvl > g_log_level) {
        return;
    }
    va_list args;
    fprintf(stderr, "[%s] ", tag);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fputc('\n', stderr);
}

#define LOG_ERR(...)  uvr_log(LOG_ERR_L,  "ERROR", __VA_ARGS__)
#define LOG_WARN(...) uvr_log(LOG_WARN_L, "WARN ", __VA_ARGS__)
#define LOG_INFO(...) uvr_log(LOG_INFO_L, "INFO ", __VA_ARGS__)
#define LOG_DBG(...)  uvr_log(LOG_DBG_L,  "DEBUG", __VA_ARGS__)

typedef enum {
    RC_MOTORWAY    = 0,
    RC_TRUNK       = 1,
    RC_PRIMARY     = 2,
    RC_SECONDARY   = 3,
    RC_TERTIARY    = 4,
    RC_RESIDENTIAL = 5,
    RC_UNCLASSIFIED= 6
} RoadClass;

typedef enum {
    NODE_NONE      = 0,
    NODE_DEAD_TREE = (1u << 0), /* Peeled 1-core cul-de-sac node      */
    NODE_JUNCTION  = (1u << 1), /* Degree != 2 intersection survivor  */
    NODE_HAS_TURN_RES = (1u << 2)
} NodeFlags;

typedef enum {
    EDGE_NONE       = 0,
    EDGE_TOLL       = (1u << 0),
    EDGE_CONTRACTED = (1u << 1), /* Collapsed macro-vector segment     */
    EDGE_BRIDGE     = (1u << 2),
    EDGE_DEAD_TARGET = (1u << 4), /* edge.to is a peeled dead-tree node */
    EDGE_DEAD_SOURCE = (1u << 5)  /* edge.from is a peeled dead-tree node */
} EdgeFlags;

typedef enum {
    MODE_VEC_BIDIR   = 0, /* Bidirectional 8D Vector-Corridor Search   */
    MODE_VEC_UNIDIR  = 1, /* Unidirectional 8D Vector-Corridor Search  */
    MODE_DIJKSTRA    = 2, /* True uninformed Dijkstra (travel time)    */
    MODE_ASTAR       = 3, /* ALT-guided A* on travel time (efficient)  */
    MODE_CCH         = 4  /* Contraction Hierarchy (CCH) query         */
} SearchMode;

/* 5D Preference Query Vector paired with 3D runtime spatial context -> 8D */
typedef struct {
    float w_align;    /* Spatial dot-product + corridor alignment weight */
    float w_speed;    /* Speed capability feature weight                 */
    float w_class;    /* Road hierarchy (motorway/trunk/primary) weight  */
    float w_straight; /* Macro-vector straightness (tortuosity) weight   */
    float w_smooth;   /* Smoothness / toll-free preference weight        */
    int   avoid_tolls;
    int   pref_id;    /* 0..4 built-in profile, -1 custom weights */
} VmrPreference;

typedef struct __attribute__((aligned(32))) {
    uint32_t to[4];
    float ux[4], uy[4];
    float dx[4], dy[4];
    float vx[4], vy[4];
    float time[4];
    float attr[5][4];   /* per built-in profile, lane-contiguous */
    uint8_t flags[4];
    uint8_t pad[12];
} EdgeBlock4;

typedef struct {
    SearchMode    mode;
    VmrPreference pref;
    int           allow_cross_track_prune;
} VmrQuery;

/* ========================================================================= */
/* SECTION 2: TRACKED ALLOCATOR & RAW POINTER UTILITIES                      */
/* ========================================================================= */

typedef struct {
    size_t magic;
    size_t bytes;
} MemHdr;

#define MEM_MAGIC ((size_t)0x555652324D454D31ULL)

static struct {
    size_t allocs;
    size_t frees;
    size_t active;
    size_t cur_bytes;
    size_t peak_bytes;
} g_mem = { 0, 0, 0, 0, 0 };

static void *xmalloc(size_t sz) {
    if (sz == 0) {
        sz = 1;
    }
    MemHdr *h = (MemHdr *)malloc(sizeof(MemHdr) + sz);
    if (h == NULL) {
        LOG_ERR("xmalloc failed for %zu bytes", sz);
        return NULL;
    }
    h->magic = MEM_MAGIC;
    h->bytes = sz;
    g_mem.allocs++;
    g_mem.active++;
    g_mem.cur_bytes += sz;
    if (g_mem.cur_bytes > g_mem.peak_bytes) {
        g_mem.peak_bytes = g_mem.cur_bytes;
    }
    return (void *)(h + 1);
}

static void *xcalloc(size_t count, size_t elem) {
    if (count != 0 && elem > SIZE_MAX / count) {
        return NULL;
    }
    size_t total = count * elem;
    void *p = xmalloc(total);
    if (p != NULL) {
        memset(p, 0, total);
    }
    return p;
}

static void *xrealloc(void *ptr, size_t new_sz) {
    if (ptr == NULL) {
        return xmalloc(new_sz);
    }
    if (new_sz == 0) {
        new_sz = 1;
    }
    MemHdr *old_h = ((MemHdr *)ptr) - 1;
    if (old_h->magic != MEM_MAGIC) {
        LOG_ERR("Corrupted header in xrealloc");
        abort();
    }
    size_t old_sz = old_h->bytes;
    MemHdr *new_h = (MemHdr *)realloc(old_h, sizeof(MemHdr) + new_sz);
    if (new_h == NULL) {
        LOG_ERR("xrealloc failed for %zu bytes", new_sz);
        return NULL;
    }
    new_h->magic = MEM_MAGIC;
    new_h->bytes = new_sz;
    g_mem.cur_bytes = (g_mem.cur_bytes - old_sz) + new_sz;
    if (g_mem.cur_bytes > g_mem.peak_bytes) {
        g_mem.peak_bytes = g_mem.cur_bytes;
    }
    return (void *)(new_h + 1);
}

static void xfree(void *ptr) {
    if (ptr == NULL) {
        return;
    }
    MemHdr *h = ((MemHdr *)ptr) - 1;
    if (h->magic != MEM_MAGIC) {
        LOG_ERR("Double or invalid free detected");
        abort();
    }
    h->magic = 0;
    g_mem.frees++;
    if (g_mem.active > 0) {
        g_mem.active--;
    }
    if (g_mem.cur_bytes >= h->bytes) {
        g_mem.cur_bytes -= h->bytes;
    }
    free(h);
}

static void *xaligned_alloc(size_t align, size_t size) {
    if (align == 0 || (align & (align - 1U)) != 0) {
        return NULL;
    }
    void *raw = xmalloc(size + align + sizeof(void *));
    if (!raw) {
        return NULL;
    }
    uintptr_t p = ((uintptr_t)raw + sizeof(void *) + (uintptr_t)align - 1U) &
                  ~((uintptr_t)align - 1U);
    ((void **)p)[-1] = raw;
    return (void *)p;
}

static void xfree_aligned(void *ptr) {
    if (!ptr) {
        return;
    }
    xfree(((void **)ptr)[-1]);
}

static int c_stricmp(const char *a, const char *b) {
    while (*a != '\0' && *b != '\0') {
        int ca = tolower((unsigned char)*a);
        int cb = tolower((unsigned char)*b);
        if (ca != cb) {
            return ca - cb;
        }
        a++;
        b++;
    }
    return tolower((unsigned char)*a) - tolower((unsigned char)*b);
}

/* ========================================================================= */
/* SECTION 3: UK CONFORMAL PROJECTION, HASH MAP & MIN-HEAP                   */
/* ========================================================================= */

static void uk_project_tm(double lat_deg, double lon_deg, float *out_x, float *out_y) {
    double phi  = lat_deg * DEG_TO_RAD;
    double phi0 = UK_TM_LAT0_DEG * DEG_TO_RAD;
    double dlam = (lon_deg - UK_TM_LON0_DEG) * DEG_TO_RAD;

    double b = cos(phi) * sin(dlam);
    if (b > 0.999999) {
        b = 0.999999;
    } else if (b < -0.999999) {
        b = -0.999999;
    }
    double x = UK_TM_K0 * EARTH_RADIUS_M * atanh(b);
    double y = UK_TM_K0 * EARTH_RADIUS_M * (atan2(tan(phi), cos(dlam)) - phi0);

    *out_x = (float)x;
    *out_y = (float)y;
}

static double haversine_m(double lat1, double lon1, double lat2, double lon2) {
    double p1 = lat1 * DEG_TO_RAD;
    double p2 = lat2 * DEG_TO_RAD;
    double dp = (lat2 - lat1) * DEG_TO_RAD;
    double dl = (lon2 - lon1) * DEG_TO_RAD;
    double sd = sin(dp * 0.5);
    double sl = sin(dl * 0.5);
    double a  = (sd * sd) + (cos(p1) * cos(p2) * sl * sl);
    if (a > 1.0) {
        a = 1.0;
    }
    return 2.0 * EARTH_RADIUS_M * atan2(sqrt(a), sqrt(1.0 - a));
}

/* Compact 16-byte open-addressing hash map (int64_t OSM ID -> uint32_t idx) */
typedef struct {
    int64_t  key;
    uint32_t val;
    uint32_t used;
} IdMapSlot;

typedef struct {
    IdMapSlot *slots;
    size_t     cap;
    size_t     count;
} IdMap;

static uint64_t mix64(int64_t k) {
    uint64_t x = (uint64_t)k;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static int idmap_init(IdMap *m, size_t min_cap) {
    size_t cap = 64;
    while (cap < min_cap * 2) {
        cap <<= 1;
    }
    m->slots = (IdMapSlot *)xcalloc(cap, sizeof(IdMapSlot));
    if (m->slots == NULL) {
        m->cap = 0;
        m->count = 0;
        return -1;
    }
    m->cap = cap;
    m->count = 0;
    return 0;
}

static void idmap_free(IdMap *m) {
    xfree(m->slots);
    m->slots = NULL;
    m->cap = 0;
    m->count = 0;
}

static int idmap_rehash(IdMap *m, size_t new_cap) {
    IdMapSlot *ns = (IdMapSlot *)xcalloc(new_cap, sizeof(IdMapSlot));
    if (ns == NULL) {
        return -1;
    }
    size_t mask = new_cap - 1;
    IdMapSlot *p = m->slots;
    IdMapSlot *end = m->slots + m->cap;
    while (p < end) {
        if (p->used) {
            size_t i = (size_t)(mix64(p->key) & mask);
            IdMapSlot *t = ns + i;
            while (t->used) {
                i = (i + 1) & mask;
                t = ns + i;
            }
            *t = *p;
        }
        p++;
    }
    xfree(m->slots);
    m->slots = ns;
    m->cap = new_cap;
    return 0;
}

static int idmap_put(IdMap *m, int64_t key, uint32_t val) {
    if (m->cap == 0) {
        if (idmap_init(m, 1024) != 0) {
            return -1;
        }
    }
    if ((m->count + 1) * 10 >= m->cap * 7) {
        if (idmap_rehash(m, m->cap * 2) != 0) {
            return -1;
        }
    }
    size_t mask = m->cap - 1;
    size_t i = (size_t)(mix64(key) & mask);
    IdMapSlot *s = m->slots + i;
    while (s->used) {
        if (s->key == key) {
            s->val = val;
            return 0;
        }
        i = (i + 1) & mask;
        s = m->slots + i;
    }
    s->key = key;
    s->val = val;
    s->used = 1;
    m->count++;
    return 0;
}

static uint32_t idmap_get(const IdMap *m, int64_t key) {
    if (m->cap == 0 || m->slots == NULL) {
        return VMR_INVALID;
    }
    size_t mask = m->cap - 1;
    size_t i = (size_t)(mix64(key) & mask);
    const IdMapSlot *s = m->slots + i;
    size_t probes = 0;
    while (s->used && probes < m->cap) {
        if (s->key == key) {
            return s->val;
        }
        i = (i + 1) & mask;
        s = m->slots + i;
        probes++;
    }
    return VMR_INVALID;
}

/* MinHeap using raw pointer arithmetic */
typedef struct {
    uint32_t node;
    float    prio;
} HeapEntry;

typedef struct {
    HeapEntry *buf;
    size_t     size;
    size_t     cap;
} MinHeap;

static int heap_init(MinHeap *h, size_t cap) {
    if (cap < 64) {
        cap = 64;
    }
    h->buf = (HeapEntry *)xmalloc(cap * sizeof(HeapEntry));
    if (h->buf == NULL) {
        h->size = 0;
        h->cap = 0;
        return -1;
    }
    h->size = 0;
    h->cap = cap;
    return 0;
}

static void heap_free(MinHeap *h) {
    xfree(h->buf);
    h->buf = NULL;
    h->size = 0;
    h->cap = 0;
}

/* Drop stale (superceded) entries, keeping only the minimum-priority entry per
 * node, then re-heapify.  Prevents unbounded growth from lazy deletion. */
static void heap_compact_min(MinHeap *h, uint32_t n) {
    size_t sz = h->size;
    if (sz <= (size_t)n * 2U || sz == 0) {
        return;
    }
    float *bp = (float *)xmalloc((size_t)n * sizeof(float));
    uint8_t *seen = (uint8_t *)xcalloc(n, 1);
    if (!bp || !seen) {
        xfree(bp); xfree(seen);
        return;
    }
    for (uint32_t i = 0; i < n; i++) {
        bp[i] = FLT_MAX;
    }
    HeapEntry *b = h->buf;
    for (size_t i = 0; i < sz; i++) {
        uint32_t nd = b[i].node;
        if (nd < n && b[i].prio < bp[nd]) {
            bp[nd] = b[i].prio;
        }
    }
    size_t w = 0;
    for (size_t i = 0; i < sz; i++) {
        uint32_t nd = b[i].node;
        if (nd < n && !seen[nd] && b[i].prio == bp[nd]) {
            seen[nd] = 1;
            b[w++] = b[i];
        }
    }
    h->size = w;
    xfree(bp);
    xfree(seen);
    /* re-heapify (bottom-up) */
    for (size_t i = h->size / 2; i-- > 0;) {
        size_t idx = i;
        while (1) {
            size_t child = idx * 4 + 1;
            if (child >= h->size) {
                break;
            }
            size_t end = child + 4;
            if (end > h->size) end = h->size;
            size_t best = child;
            for (size_t c = child + 1; c < end; c++) {
                if (b[c].prio < b[best].prio) best = c;
            }
            if (b[best].prio < b[idx].prio) {
                HeapEntry tmp = b[idx]; b[idx] = b[best]; b[best] = tmp;
                idx = best;
            } else {
                break;
            }
        }
    }
}

static int heap_push(MinHeap *h, uint32_t node, float prio) {
    if (h->size >= h->cap) {
        size_t ncap = h->cap * 2;
        HeapEntry *nb = (HeapEntry *)xrealloc(h->buf, ncap * sizeof(HeapEntry));
        if (nb == NULL) {
            return -1;
        }
        h->buf = nb;
        h->cap = ncap;
    }
    size_t i = h->size++;
    HeapEntry *base = h->buf;
    *(base + i) = (HeapEntry){ node, prio };

    while (i > 0) {
        size_t p = (i - 1) >> 2;
        HeapEntry *ci = base + i;
        HeapEntry *cp = base + p;
        if (ci->prio < cp->prio) {
            HeapEntry tmp = *ci;
            *ci = *cp;
            *cp = tmp;
            i = p;
        } else {
            break;
        }
    }
    return 0;
}

static int heap_pop(MinHeap *h, HeapEntry *out) {
    if (h->size == 0) {
        return 0;
    }
    HeapEntry *base = h->buf;
    *out = *base;
    h->size--;
    if (h->size > 0) {
        *base = *(base + h->size);
        size_t i = 0;
        size_t n = h->size;
        while (1) {
            size_t child = i * 4 + 1;
            if (child >= n) break;
            size_t end = child + 4;
            if (end > n) end = n;
            size_t best = child;
            for (size_t c = child + 1; c < end; c++) {
                if ((base + c)->prio < (base + best)->prio) best = c;
            }
            if ((base + best)->prio < (base + i)->prio) {
                HeapEntry tmp = *(base + i);
                *(base + i) = *(base + best);
                *(base + best) = tmp;
                i = best;
            } else {
                break;
            }
        }
    }
    return 1;
}

/* Monotone bucket queue (Dial): A* keys are non-decreasing with the consistent
 * ALT heuristic, so integer buckets (1 s) replace the 4-ary heap. */
#define BUCKET_NB (1u << 20)
typedef struct {
    uint32_t *head;   /* [BUCKET_NB] node or VMR_INVALID */
    uint32_t *next;   /* [n] next node in the same bucket */
    uint32_t *prev;   /* [n] prev node in the same bucket */
    uint32_t *bucket; /* [n] current bucket or VMR_INVALID */
    uint32_t cur;
    uint32_t maxb;
    uint32_t nb;
    uint32_t size;
    size_t cap;
} BucketQueue;

static inline void bq_clear(BucketQueue *q);

static __attribute__((unused)) int bq_init(BucketQueue *q, size_t n) {
    q->head = (uint32_t *)xmalloc((size_t)BUCKET_NB * sizeof(uint32_t));
    q->next = (uint32_t *)xmalloc((n ? n : 1U) * sizeof(uint32_t));
    q->prev = (uint32_t *)xmalloc((n ? n : 1U) * sizeof(uint32_t));
    q->bucket = (uint32_t *)xmalloc((n ? n : 1U) * sizeof(uint32_t));
    if (!q->head || !q->next || !q->prev || !q->bucket) {
        xfree(q->head);
        xfree(q->next);
        xfree(q->prev);
        xfree(q->bucket);
        q->head = NULL;
        q->next = NULL;
        q->prev = NULL;
        q->bucket = NULL;
        return -1;
    }
    q->nb = BUCKET_NB;
    q->cap = n ? n : 1U;
    bq_clear(q);
    return 0;
}

static __attribute__((unused)) void bq_free(BucketQueue *q) {
    xfree(q->head);
    xfree(q->next);
    xfree(q->prev);
    xfree(q->bucket);
    q->head = NULL;
    q->next = NULL;
    q->prev = NULL;
    q->bucket = NULL;
}

static inline void bq_clear(BucketQueue *q) {
    memset(q->head, 0xFF, (size_t)q->nb * sizeof(uint32_t));
    memset(q->bucket, 0xFF, (size_t)q->cap * sizeof(uint32_t));
    q->cur = 0;
    q->maxb = 0;
    q->size = 0;
}

static inline void bq_remove(BucketQueue *q, uint32_t node) {
    uint32_t b = q->bucket[node];
    if (b == VMR_INVALID) {
        return;
    }
    uint32_t p = q->prev[node], n = q->next[node];
    if (p != VMR_INVALID) {
        q->next[p] = n;
    } else {
        q->head[b] = n;
    }
    if (n != VMR_INVALID) {
        q->prev[n] = p;
    }
    q->bucket[node] = VMR_INVALID;
    q->size--;
}

static inline int bq_push(BucketQueue *q, uint32_t node, float prio) {
    if (q->bucket[node] != VMR_INVALID) {
        bq_remove(q, node);
    }
    uint32_t b = (uint32_t)prio;
    if (b >= q->nb) {
        b = q->nb - 1U;
    }
    q->next[node] = q->head[b];
    q->prev[node] = VMR_INVALID;
    if (q->head[b] != VMR_INVALID) {
        q->prev[q->head[b]] = node;
    }
    q->head[b] = node;
    q->bucket[node] = b;
    if (b > q->maxb) {
        q->maxb = b;
    }
    q->size++;
    return 0;
}

static inline int bq_pop(BucketQueue *q, HeapEntry *out) {
    while (q->cur <= q->maxb && q->head[q->cur] == VMR_INVALID) {
        q->cur++;
    }
    if (q->cur > q->maxb) {
        q->size = 0;
        return 0;
    }
    uint32_t node = q->head[q->cur];
    q->head[q->cur] = q->next[node];
    if (q->head[q->cur] != VMR_INVALID) {
        q->prev[q->head[q->cur]] = VMR_INVALID;
    }
    q->bucket[node] = VMR_INVALID;
    q->size--;
    out->node = node;
    out->prio = (float)q->cur;
    return 1;
}

static inline float bq_top(const BucketQueue *q) {
    uint32_t b = q->cur;
    while (b <= q->maxb && q->head[b] == VMR_INVALID) {
        b++;
    }
    return (b <= q->maxb) ? (float)b : 1.0e30f;
}

/* Radix heap for monotone float A* keys (consistent ALT heuristic). */
#define RADIX_B 33
typedef struct {
    HeapEntry *a[RADIX_B];
    size_t n[RADIX_B], cap[RADIX_B];
    size_t size;
    float last;
} RadixHeap;

static __attribute__((unused)) int rh_init(RadixHeap *q) {
    memset(q, 0, sizeof(*q));
    return 0;
}

static __attribute__((unused)) void rh_free(RadixHeap *q) {
    for (int b = 0; b < RADIX_B; b++) {
        xfree(q->a[b]);
    }
    memset(q, 0, sizeof(*q));
}

static inline void rh_clear(RadixHeap *q) {
    for (int b = 0; b < RADIX_B; b++) {
        q->n[b] = 0;
    }
    q->size = 0;
    q->last = 0.0f;
}

static inline int rh_push(RadixHeap *q, uint32_t node, float prio) {
    if (q->size == 0) {
        q->last = prio;
    }
    uint32_t kb, kp;
    memcpy(&kb, &prio, 4);
    memcpy(&kp, &q->last, 4);
    uint32_t x = kb ^ kp;
    int b = (x == 0) ? 0 : (31 - __builtin_clz(x));
    if (b >= RADIX_B) {
        b = RADIX_B - 1;
    }
    if (q->n[b] >= q->cap[b]) {
        size_t nc = q->cap[b] ? q->cap[b] * 2 : 64;
        HeapEntry *na = (HeapEntry *)xrealloc(q->a[b], nc * sizeof(HeapEntry));
        if (!na) {
            return -1;
        }
        q->a[b] = na;
        q->cap[b] = nc;
    }
    q->a[b][q->n[b]++] = (HeapEntry){ node, prio };
    q->size++;
    return 0;
}

static inline int rh_pop(RadixHeap *q, HeapEntry *out) {
    if (q->size == 0) {
        return 0;
    }
    int b = 0;
    while (b < RADIX_B && q->n[b] == 0) {
        b++;
    }
    if (b == RADIX_B) {
        return 0;
    }
    size_t mi = 0;
    for (size_t i = 1; i < q->n[b]; i++) {
        if (q->a[b][i].prio < q->a[b][mi].prio) {
            mi = i;
        }
    }
    *out = q->a[b][mi];
    q->a[b][mi] = q->a[b][--q->n[b]];
    q->size--;
    q->last = out->prio;
    return 1;
}

static inline float rh_top(const RadixHeap *q) {
    return q->last;
}

/* ========================================================================= */
/* SECTION 4: CONTIGUOUS SoA GRAPH & PREPROCESSED GRAPH DATA STRUCTURES      */
/* ========================================================================= */

typedef struct {
    uint32_t from;
    uint32_t to;
    float    length_m;
    float    travel_time_s;
    float    f_speed;
    float    f_class;
    float    f_smooth;
    uint8_t  road_class;
    uint8_t  flags;
    uint8_t  alive;
} RawEdge;

typedef struct {
    int64_t *node_osm_id;
    double  *node_lat;
    double  *node_lon;
    float   *node_x;
    float   *node_y;
    uint32_t node_count;
    uint32_t node_cap;

    RawEdge *edges;
    uint32_t edge_count;
    uint32_t edge_cap;

    uint32_t *turn_via;
    uint32_t *turn_from;
    uint32_t *turn_to;
    uint8_t  *turn_mode;
    uint32_t turn_count;

    IdMap    osm_map;
} VmrRawGraph;

typedef struct __attribute__((aligned(64))) {
    float fwd[32];
    float bwd[32];
} AltNode32;

/*
 * Final Preprocessed Graph in Structure-of-Arrays (SoA) Contiguous Memory:
 *
 * Every directed edge `e` carries an 8D vector:
 *   [ dx, dy, ux, uy, f_speed, f_class, f_straight, f_smooth ]
 */
typedef struct {
    uint32_t node_count;
    uint32_t edge_count;
    uint32_t raw_node_count;
    uint32_t raw_edge_count;
    uint32_t contracted_nodes;
    uint32_t peeled_dead_nodes;

    /* Node SoA */
    int64_t  *node_osm_id;
    double   *node_lat;
    double   *node_lon;
    float    *node_x;
    float    *node_y;
    uint32_t *node_comp_id;
    uint8_t  *node_flags;

    /* Forward CSR offsets (len = node_count + 1) */
    uint32_t *fwd_offsets;

    /* Forward Edge 8D Vector SoA (len = edge_count) */
    uint32_t *edge_from;
    uint32_t *edge_to;
    float    *edge_dx;
    float    *edge_dy;
    float    *edge_ux;         /* Precomputed unit dx */
    float    *edge_uy;         /* Precomputed unit dy */
    float    *edge_f_speed;    /* [0..1] normalized speed capability */
    float    *edge_f_class;    /* [0..1] highway hierarchy           */
    float    *edge_f_straight; /* (0..1] macro-vector straightness   */
    float    *edge_f_smooth;   /* [0..1] access/toll smoothness      */
    float    *edge_length_m;
    float    *edge_time_s;
    uint8_t  *edge_flags;

    /* Reverse CSR for bidirectional search (len = node_count + 1 / edge_count) */
    uint32_t *rev_offsets;
    uint32_t *rev_src_node;
    uint32_t *rev_fwd_eid;

    /* 128x128 Spatial Snap Grid (maps lat/lon -> nearest intersection in ~O(1)) */
    float     grid_min_x;
    float     grid_max_x;
    float     grid_min_y;
    float     grid_max_y;
    uint32_t *grid_offsets;    /* len = SPATIAL_GRID_DIM * SPATIAL_GRID_DIM + 1 */
    uint32_t *grid_nodes;      /* len = node_count */

    /* ALT landmarks: L sources, dist(landmark -> node) and dist(node -> landmark) */
    uint32_t  landmark_count;
    uint32_t *landmark_node;
    AltNode32 *alt_nodes;

    /* Turn restrictions: per-via forbidden (from,to) pairs */
    uint32_t *turn_offsets;    /* [node_count + 1] */
    uint32_t *turn_from;
    uint32_t *turn_to;
    uint8_t  *turn_mode;       /* 0=forbid */
    uint32_t  turn_count;

    /* Vector-aware contraction hierarchy (shortcuts carry the same 8D features) */
    uint32_t  ch_built;
    uint32_t  cch_built;
    uint32_t *ch_level;
    uint32_t  ch_edge_count;
    uint32_t *ch_from;
    uint32_t *ch_to;
    float    *ch_time;
    uint32_t *ch_fwd_off;
    uint32_t *ch_rev_off;
    uint32_t *ch_rev_src;
    uint32_t *ch_rev_eid;
    uint32_t *ch_via;
    uint32_t *ch_child1;
    uint32_t *ch_child2;

    float *edge_attr_scale[5];
    float  attr_min[5]; /* per-profile min(edge_attr_scale) (admissible ALT scale) */

    EdgeBlock4 *rev_blocks;
    uint32_t *rev_block_off;

    /* Precomputed two-hop continuations per edge (sum vector + inverse norm) */
    uint32_t *th_off;   /* [edge_count + 1] */
    float    *th_sx;    /* [th_count] */
    float    *th_sy;
    float    *th_inv;
    uint32_t  th_count;
    uint32_t *rth_off;  /* reverse edge -> continuations */
    float    *rth_sx;
    float    *rth_sy;
    float    *rth_inv;
    uint32_t  rth_count;

    IdMap     osm_map;
    void     *map;
    HANDLE    map_h;
} VmrGraph;

typedef struct {
    int        found;
    SearchMode mode_used;
    uint32_t   start_idx;
    uint32_t   goal_idx;

    uint32_t  *nodes;          /* Sequence start -> ... -> goal */
    uint32_t  *edges;          /* Sequence of edge indices      */
    uint32_t   hop_count;

    double     total_distance_m;
    double     direct_distance_m;
    double     est_time_s;
    double     detour_ratio;
    double     mean_cos_sim;
    double     cum_vec_cos_sim;
    double     mean_8d_alignment;
    uint32_t   nodes_expanded;
    double     elapsed_ms;
} VmrResult;

static void raw_graph_init(VmrRawGraph *rg) {
    memset(rg, 0, sizeof(*rg));
}

static void raw_graph_free(VmrRawGraph *rg) {
    xfree(rg->node_osm_id);
    xfree(rg->node_lat);
    xfree(rg->node_lon);
    xfree(rg->node_x);
    xfree(rg->node_y);
    xfree(rg->edges);
    xfree(rg->turn_via);
    xfree(rg->turn_from);
    xfree(rg->turn_to);
    xfree(rg->turn_mode);
    idmap_free(&rg->osm_map);
    raw_graph_init(rg);
}

static int raw_graph_reserve_nodes(VmrRawGraph *rg, uint32_t cap) {
    if (cap <= rg->node_cap) {
        return 0;
    }
    uint32_t ncap = (rg->node_cap == 0) ? 1024U : rg->node_cap;
    while (ncap < cap) {
        ncap *= 2U;
    }
    int64_t *oids = (int64_t *)xrealloc(rg->node_osm_id, (size_t)ncap * sizeof(int64_t));
    if (!oids) {
        return -1;
    }
    rg->node_osm_id = oids;
    double *lats = (double *)xrealloc(rg->node_lat, (size_t)ncap * sizeof(double));
    if (!lats) {
        return -1;
    }
    rg->node_lat = lats;
    double *lons = (double *)xrealloc(rg->node_lon, (size_t)ncap * sizeof(double));
    if (!lons) {
        return -1;
    }
    rg->node_lon = lons;
    float *xs = (float *)xrealloc(rg->node_x, (size_t)ncap * sizeof(float));
    if (!xs) {
        return -1;
    }
    rg->node_x = xs;
    float *ys = (float *)xrealloc(rg->node_y, (size_t)ncap * sizeof(float));
    if (!ys) {
        return -1;
    }
    rg->node_y = ys;
    rg->node_cap = ncap;
    return 0;
}

static int raw_graph_add_node(VmrRawGraph *rg, int64_t osm_id, double lat, double lon) {
    if (raw_graph_reserve_nodes(rg, rg->node_count + 1U) != 0) {
        return -1;
    }
    uint32_t idx = rg->node_count;
    *(rg->node_osm_id + idx) = osm_id;
    *(rg->node_lat + idx)    = lat;
    *(rg->node_lon + idx)    = lon;
    uk_project_tm(lat, lon, rg->node_x + idx, rg->node_y + idx);
    if (idmap_put(&rg->osm_map, osm_id, idx) != 0) {
        return -1;
    }
    rg->node_count++;
    return 0;
}

static int raw_graph_add_edge(VmrRawGraph *rg, uint32_t from, uint32_t to,
                              float len_m, float time_s, float f_speed,
                              float f_class, float f_smooth,
                              uint8_t road_class, uint8_t flags) {
    if (from == to || from >= rg->node_count || to >= rg->node_count) {
        return 0;
    }
    if (rg->edge_count >= rg->edge_cap) {
        uint32_t ncap = (rg->edge_cap == 0) ? 2048U : rg->edge_cap * 2U;
        RawEdge *arr = (RawEdge *)xrealloc(rg->edges, (size_t)ncap * sizeof(RawEdge));
        if (arr == NULL) {
            return -1;
        }
        rg->edges = arr;
        rg->edge_cap = ncap;
    }
    RawEdge *e = rg->edges + rg->edge_count;
    e->from = from;
    e->to = to;
    e->length_m = (len_m > 0.2f) ? len_m : 0.2f;
    e->travel_time_s = (time_s > 0.01f) ? time_s : 0.01f;
    e->f_speed = f_speed;
    e->f_class = f_class;
    e->f_smooth = f_smooth;
    e->road_class = road_class;
    e->flags = flags;
    e->alive = 1;
    rg->edge_count++;
    return 0;
}

static void vmr_unmap_file(void *base, HANDLE h) {
    if (base) {
        UnmapViewOfFile(base);
    }
    if (h) {
        CloseHandle(h);
    }
}

static void vmr_graph_init(VmrGraph *g) {
    memset(g, 0, sizeof(*g));
    g->map = NULL;
    g->map_h = NULL;
    for (int i = 0; i < 5; i++) g->attr_min[i] = 1.0f;
}

/* Compute per-profile min(edge_attr_scale) (admissible ALT scale). */
static void vmr_compute_attr_min(VmrGraph *g) {
    for (uint32_t p = 0; p < 5; p++) {
        const float *a = g->edge_attr_scale[p];
        if (!a) { g->attr_min[p] = 1.0f; continue; }
        float mn = 1.0f;
        for (uint32_t e = 0; e < g->edge_count; e++) {
            if (a[e] < mn) mn = a[e];
        }
        g->attr_min[p] = mn;
    }
}

static void vmr_graph_free(VmrGraph *g) {
    if (!g->map) xfree(g->node_osm_id);
    if (!g->map) xfree(g->node_lat);
    if (!g->map) xfree(g->node_lon);
    if (!g->map) xfree(g->node_x);
    if (!g->map) xfree(g->node_y);
    if (!g->map) xfree(g->node_comp_id);
    if (!g->map) xfree(g->node_flags);
    if (!g->map) xfree(g->fwd_offsets);
    if (!g->map) xfree(g->edge_from);
    if (!g->map) xfree(g->edge_to);
    if (!g->map) xfree(g->edge_dx);
    if (!g->map) xfree(g->edge_dy);
    if (!g->map) xfree(g->edge_ux);
    if (!g->map) xfree(g->edge_uy);
    if (!g->map) xfree(g->edge_f_speed);
    if (!g->map) xfree(g->edge_f_class);
    if (!g->map) xfree(g->edge_f_straight);
    if (!g->map) xfree(g->edge_f_smooth);
    if (!g->map) xfree(g->edge_length_m);
    if (!g->map) xfree(g->edge_time_s);
    if (!g->map) xfree(g->edge_flags);
    if (!g->map) xfree(g->rev_offsets);
    if (!g->map) xfree(g->rev_src_node);
    if (!g->map) xfree(g->rev_fwd_eid);
    if (!g->map) xfree(g->grid_offsets);
    if (!g->map) xfree(g->grid_nodes);
    if (!g->map) xfree(g->landmark_node);
    if (!g->map) xfree_aligned(g->alt_nodes);
    xfree(g->turn_offsets);
    xfree(g->turn_from);
    xfree(g->turn_to);
    xfree(g->turn_mode);
    if (!g->map) {
        xfree(g->ch_level);
        xfree(g->ch_from);
        xfree(g->ch_to);
        xfree(g->ch_time);
        xfree(g->ch_fwd_off);
        xfree(g->ch_rev_off);
        xfree(g->ch_rev_src);
        xfree(g->ch_rev_eid);
        xfree(g->ch_via);
        xfree(g->ch_child1);
        xfree(g->ch_child2);
    }
    for (int i = 0; i < 5; i++) {
        xfree(g->edge_attr_scale[i]);
    }
    xfree_aligned(g->rev_blocks);
    xfree(g->rev_block_off);
    if (!g->map) {
        xfree(g->th_off); xfree(g->th_sx); xfree(g->th_sy); xfree(g->th_inv);
        xfree(g->rth_off); xfree(g->rth_sx); xfree(g->rth_sy); xfree(g->rth_inv);
    }
    g->th_count = 0;
    g->rth_count = 0;
    idmap_free(&g->osm_map);
    vmr_unmap_file(g->map, g->map_h);
    g->map = NULL;
    g->map_h = NULL;
    vmr_graph_init(g);
}

static void vmr_result_init(VmrResult *r) {
    memset(r, 0, sizeof(*r));
    r->start_idx = VMR_INVALID;
    r->goal_idx = VMR_INVALID;
}

static void vmr_result_free(VmrResult *r) {
    xfree(r->nodes);
    xfree(r->edges);
    vmr_result_init(r);
}

/* ========================================================================= */
/* SECTION 5: O(V + E) PREPROCESSING (CONTRACTION, PEELING, SoA, GRID)       */
/* ========================================================================= */

static int cmp_raw_edge(const void *a, const void *b) {
    const RawEdge *ea = (const RawEdge *)a;
    const RawEdge *eb = (const RawEdge *)b;
    if (ea->from < eb->from) return -1;
    if (ea->from > eb->from) return 1;
    if (ea->to < eb->to) return -1;
    if (ea->to > eb->to) return 1;
    return 0;
}

static uint32_t grid_cell_for_xy(const VmrGraph *g, float x, float y) {
    float span_x = g->grid_max_x - g->grid_min_x;
    float span_y = g->grid_max_y - g->grid_min_y;
    if (span_x < 1.0f) span_x = 1.0f;
    if (span_y < 1.0f) span_y = 1.0f;

    int gx = (int)(((x - g->grid_min_x) / span_x) * (float)(SPATIAL_GRID_DIM - 1));
    int gy = (int)(((y - g->grid_min_y) / span_y) * (float)(SPATIAL_GRID_DIM - 1));
    if (gx < 0) gx = 0; else if (gx >= SPATIAL_GRID_DIM) gx = SPATIAL_GRID_DIM - 1;
    if (gy < 0) gy = 0; else if (gy >= SPATIAL_GRID_DIM) gy = SPATIAL_GRID_DIM - 1;
    return (uint32_t)(gy * SPATIAL_GRID_DIM + gx);
}

static int vmr_build_spatial_grid(VmrGraph *g) {
    uint32_t n = g->node_count;
    if (n == 0) {
        return -1;
    }
    float min_x = FLT_MAX, max_x = -FLT_MAX;
    float min_y = FLT_MAX, max_y = -FLT_MAX;
    const float *xp = g->node_x;
    const float *yp = g->node_y;
    for (uint32_t i = 0; i < n; i++) {
        float x = *(xp + i);
        float y = *(yp + i);
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }
    g->grid_min_x = min_x;
    g->grid_max_x = max_x;
    g->grid_min_y = min_y;
    g->grid_max_y = max_y;

    uint32_t cells = (uint32_t)(SPATIAL_GRID_DIM * SPATIAL_GRID_DIM);
    g->grid_offsets = (uint32_t *)xcalloc((size_t)cells + 1U, sizeof(uint32_t));
    g->grid_nodes   = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint32_t *write_pos = (uint32_t *)xmalloc((size_t)cells * sizeof(uint32_t));
    if (!g->grid_offsets || !g->grid_nodes || !write_pos) {
        xfree(write_pos);
        return -1;
    }

    for (uint32_t i = 0; i < n; i++) {
        uint32_t c = grid_cell_for_xy(g, *(xp + i), *(yp + i));
        (*(g->grid_offsets + c + 1U))++;
    }
    for (uint32_t c = 1; c <= cells; c++) {
        *(g->grid_offsets + c) += *(g->grid_offsets + c - 1U);
    }
    memcpy(write_pos, g->grid_offsets, (size_t)cells * sizeof(uint32_t));
    for (uint32_t i = 0; i < n; i++) {
        uint32_t c = grid_cell_for_xy(g, *(xp + i), *(yp + i));
        uint32_t pos = (*(write_pos + c))++;
        *(g->grid_nodes + pos) = i;
    }
    xfree(write_pos);
    return 0;
}

static int vmr_build_reverse_csr_and_components(VmrGraph *g) {
    uint32_t n = g->node_count;
    uint32_t m = g->edge_count;

    g->rev_offsets  = (uint32_t *)xcalloc((size_t)n + 1U, sizeof(uint32_t));
    g->rev_src_node = (uint32_t *)xmalloc((size_t)(m > 0 ? m : 1) * sizeof(uint32_t));
    g->rev_fwd_eid  = (uint32_t *)xmalloc((size_t)(m > 0 ? m : 1) * sizeof(uint32_t));
    uint32_t *rpos  = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    if (!g->rev_offsets || !g->rev_src_node || !g->rev_fwd_eid || !rpos) {
        xfree(rpos);
        return -1;
    }

    const uint32_t *to_p   = g->edge_to;
    const uint32_t *from_p = g->edge_from;
    for (uint32_t e = 0; e < m; e++) {
        (*(g->rev_offsets + *(to_p + e) + 1U))++;
    }
    for (uint32_t i = 1; i <= n; i++) {
        *(g->rev_offsets + i) += *(g->rev_offsets + i - 1U);
    }
    memcpy(rpos, g->rev_offsets, (size_t)n * sizeof(uint32_t));
    for (uint32_t e = 0; e < m; e++) {
        uint32_t dst = *(to_p + e);
        uint32_t idx = (*(rpos + dst))++;
        *(g->rev_src_node + idx) = *(from_p + e);
        *(g->rev_fwd_eid + idx)  = e;
    }
    xfree(rpos);

    /* Connected components via undirected BFS (fwd + rev) for O(1) island checks */
    g->node_comp_id = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint32_t *queue = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    if (!g->node_comp_id || !queue) {
        xfree(queue);
        return -1;
    }
    for (uint32_t i = 0; i < n; i++) {
        *(g->node_comp_id + i) = VMR_INVALID;
    }

    uint32_t comp = 0;
    for (uint32_t s = 0; s < n; s++) {
        if (*(g->node_comp_id + s) != VMR_INVALID) {
            continue;
        }
        uint32_t head = 0, tail = 0;
        *(g->node_comp_id + s) = comp;
        *(queue + tail++) = s;
        while (head < tail) {
            uint32_t u = *(queue + head++);
            uint32_t fb = *(g->fwd_offsets + u);
            uint32_t fe = *(g->fwd_offsets + u + 1U);
            for (uint32_t ei = fb; ei < fe; ei++) {
                uint32_t v = *(to_p + ei);
                if (*(g->node_comp_id + v) == VMR_INVALID) {
                    *(g->node_comp_id + v) = comp;
                    *(queue + tail++) = v;
                }
            }
            uint32_t rb = *(g->rev_offsets + u);
            uint32_t re = *(g->rev_offsets + u + 1U);
            for (uint32_t ri = rb; ri < re; ri++) {
                uint32_t v = *(g->rev_src_node + ri);
                if (*(g->node_comp_id + v) == VMR_INVALID) {
                    *(g->node_comp_id + v) = comp;
                    *(queue + tail++) = v;
                }
            }
        }
        comp++;
    }
    xfree(queue);
    LOG_INFO("Connected components labeled: %u component(s) (mainland/islands)", comp);
    return 0;
}

/* ------------------------------------------------------------------------- */
/* Vector-aware contraction hierarchy (CH)                                    */
/* ------------------------------------------------------------------------- */

typedef struct {
    uint32_t from, to;
    float length_m, time_s;
    float f_speed, f_class, f_smooth;
    float dx, dy;
    uint8_t road_class, flags, alive;
    uint32_t via;
    uint32_t child1;
    uint32_t child2;
} ChEdge;

typedef struct {
    ChEdge *edges;
    size_t ecount, ecap;
    uint32_t *out_head, *in_head;
    uint32_t *out_next, *in_next;
    uint32_t *out_deg, *in_deg;
    uint32_t n;
} ChAdj;

typedef struct {
    uint64_t *keys;
    uint8_t *used;
    uint32_t *vals;
    size_t cap, count;
} ChPairSet;

typedef struct {
    uint32_t from, to;
    float length_m, time_s;
    float f_speed, f_class, f_smooth;
    uint8_t flags;
    uint32_t via;
    uint32_t child1;
    uint32_t child2;
} ChPacked;

static int ch_adj_init(ChAdj *a, uint32_t n, uint32_t initial) {
    memset(a, 0, sizeof(*a));
    a->n = n;
    size_t cap = initial ? initial : 16;
    a->edges = (ChEdge *)xmalloc(cap * sizeof(ChEdge));
    a->out_head = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    a->in_head  = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    a->out_next = (uint32_t *)xmalloc(cap * sizeof(uint32_t));
    a->in_next  = (uint32_t *)xmalloc(cap * sizeof(uint32_t));
    a->out_deg  = (uint32_t *)xcalloc(n, sizeof(uint32_t));
    a->in_deg   = (uint32_t *)xcalloc(n, sizeof(uint32_t));
    if (!a->edges || !a->out_head || !a->in_head || !a->out_next ||
        !a->in_next || !a->out_deg || !a->in_deg) {
        return -1;
    }
    a->ecap = cap;
    for (uint32_t i = 0; i < n; i++) {
        *(a->out_head + i) = VMR_INVALID;
        *(a->in_head + i)  = VMR_INVALID;
    }
    return 0;
}

static void ch_adj_free(ChAdj *a) {
    xfree(a->edges); xfree(a->out_head); xfree(a->in_head);
    xfree(a->out_next); xfree(a->in_next);
    xfree(a->out_deg); xfree(a->in_deg);
    memset(a, 0, sizeof(*a));
}

static uint32_t ch_adj_add(ChAdj *a, uint32_t from, uint32_t to,
                           float length_m, float time_s,
                           float f_speed, float f_class, float f_smooth,
                           float dx, float dy, uint8_t road_class,
                           uint8_t flags, uint32_t via,
                           uint32_t child1, uint32_t child2) {
    if (a->ecount == a->ecap) {
        size_t ncap = a->ecap * 2U;
        ChEdge *ne = (ChEdge *)xrealloc(a->edges, ncap * sizeof(ChEdge));
        if (!ne) return VMR_INVALID;
        a->edges = ne;
        uint32_t *on = (uint32_t *)xrealloc(a->out_next, ncap * sizeof(uint32_t));
        if (!on) return VMR_INVALID;
        a->out_next = on;
        uint32_t *inn = (uint32_t *)xrealloc(a->in_next, ncap * sizeof(uint32_t));
        if (!inn) return VMR_INVALID;
        a->in_next = inn;
        a->ecap = ncap;
    }
    uint32_t id = (uint32_t)a->ecount++;
    ChEdge *e = a->edges + id;
    e->from = from; e->to = to;
    e->length_m = length_m; e->time_s = time_s;
    e->f_speed = f_speed; e->f_class = f_class; e->f_smooth = f_smooth;
    e->dx = dx; e->dy = dy;
    e->road_class = road_class; e->flags = flags; e->alive = 1;
    e->via = via;
    e->child1 = child1;
    e->child2 = child2;
    *(a->out_next + id) = *(a->out_head + from);
    *(a->out_head + from) = id;
    (*(a->out_deg + from))++;
    *(a->in_next + id) = *(a->in_head + to);
    *(a->in_head + to) = id;
    (*(a->in_deg + to))++;
    return id;
}

static void ch_adj_kill(ChAdj *a, uint32_t id) {
    ChEdge *e = a->edges + id;
    if (e->alive) {
        e->alive = 0;
        (*(a->out_deg + e->from))--;
        (*(a->in_deg + e->to))--;
    }
}

static int ch_pair_init(ChPairSet *s, size_t cap) {
    memset(s, 0, sizeof(*s));
    s->cap = cap;
    s->keys = (uint64_t *)xmalloc(cap * sizeof(uint64_t));
    s->used = (uint8_t *)xcalloc(cap, 1);
    s->vals = (uint32_t *)xmalloc(cap * sizeof(uint32_t));
    return (s->keys && s->used && s->vals) ? 0 : -1;
}

static void ch_pair_free(ChPairSet *s) {
    xfree(s->keys); xfree(s->used); xfree(s->vals);
    memset(s, 0, sizeof(*s));
}

static int ch_pair_add(ChPairSet *s, uint32_t u, uint32_t v) {
    if ((s->count + 1U) * 2U >= s->cap) {
        size_t oldcap = s->cap;
        size_t ncap = oldcap * 2U;
        uint64_t *nk = (uint64_t *)xmalloc(ncap * sizeof(uint64_t));
        uint8_t *nu = (uint8_t *)xcalloc(ncap, 1);
        uint32_t *nv = (uint32_t *)xmalloc(ncap * sizeof(uint32_t));
        if (!nk || !nu || !nv) {
            xfree(nk); xfree(nu); xfree(nv);
            return -1;
        }
        size_t nmask = ncap - 1U;
        for (size_t i = 0; i < oldcap; i++) {
            if (!*(s->used + i)) continue;
            uint64_t key = *(s->keys + i);
            size_t j = (size_t)(mix64((int64_t)key) & nmask);
            while (*(nu + j)) {
                j = (j + 1U) & nmask;
            }
            *(nu + j) = 1;
            *(nk + j) = key;
            *(nv + j) = *(s->vals + i);
        }
        xfree(s->keys); xfree(s->used); xfree(s->vals);
        s->keys = nk; s->used = nu; s->vals = nv;
        s->cap = ncap;
    }
    uint64_t key = ((uint64_t)u << 32) | (uint64_t)v;
    size_t mask = s->cap - 1U;
    size_t i = (size_t)(mix64((int64_t)key) & mask);
    while (*(s->used + i)) {
        if (*(s->keys + i) == key) {
            return 0;
        }
        i = (i + 1U) & mask;
    }
    *(s->used + i) = 1;
    *(s->keys + i) = key;
    *(s->vals + i) = VMR_INVALID;
    s->count++;
    return 1;
}

static int ch_pair_get(const ChPairSet *s, uint32_t u, uint32_t v, uint32_t *eid) {
    if (s->cap == 0) return 0;
    uint64_t key = ((uint64_t)u << 32) | (uint64_t)v;
    size_t mask = s->cap - 1U;
    size_t i = (size_t)(mix64((int64_t)key) & mask);
    while (*(s->used + i)) {
        if (*(s->keys + i) == key) {
            if (eid) *eid = *(s->vals + i);
            return 1;
        }
        i = (i + 1U) & mask;
    }
    return 0;
}

static void ch_pair_set(ChPairSet *s, uint32_t u, uint32_t v, uint32_t eid) {
    uint64_t key = ((uint64_t)u << 32) | (uint64_t)v;
    size_t mask = s->cap - 1U;
    size_t i = (size_t)(mix64((int64_t)key) & mask);
    while (*(s->used + i)) {
        if (*(s->keys + i) == key) {
            *(s->vals + i) = eid;
            return;
        }
        i = (i + 1U) & mask;
    }
}

static int ch_contract_node(ChAdj *a, uint32_t u, uint8_t *contracted,
                            ChPairSet *pairs, uint32_t *level,
                            uint32_t *counter, uint32_t *cur, MinHeap *pq) {
    if (u >= a->n || *(contracted + u)) {
        return 0;
    }
    for (uint32_t e1 = *(a->in_head + u); e1 != VMR_INVALID; e1 = *(a->in_next + e1)) {
        if (!(*(a->edges + e1)).alive) continue;
        uint32_t uu = (*(a->edges + e1)).from;
        for (uint32_t e2 = *(a->out_head + u); e2 != VMR_INVALID; e2 = *(a->out_next + e2)) {
            if (!(*(a->edges + e2)).alive) continue;
            uint32_t vv = (*(a->edges + e2)).to;
            if (uu == vv || *(contracted + uu) || *(contracted + vv)) continue;
            int added = ch_pair_add(pairs, uu, vv);
            if (added < 0) return -1;
            ChEdge *x = a->edges + e1;
            ChEdge *y = a->edges + e2;
            float len = x->length_m + y->length_m;
            float t = x->time_s + y->time_s;
            float fs = (x->f_speed < y->f_speed) ? x->f_speed : y->f_speed;
            float fc = (x->f_class < y->f_class) ? x->f_class : y->f_class;
            float fm = (x->f_smooth < y->f_smooth) ? x->f_smooth : y->f_smooth;
            float dx = x->dx + y->dx;
            float dy = x->dy + y->dy;
            if (added == 0) {
                uint32_t eid;
                if (ch_pair_get(pairs, uu, vv, &eid) && t < (*(a->edges + eid)).time_s) {
                    ChEdge *ex = a->edges + eid;
                    ex->time_s = t;
                    ex->length_m = len;
                    ex->f_speed = fs;
                    ex->f_class = fc;
                    ex->f_smooth = fm;
                    ex->dx = dx;
                    ex->dy = dy;
                    ex->flags = (uint8_t)(x->flags | y->flags);
                    /* Rebind the shortcut to the shorter path's children (e1,e2
                     * are strict sub-edges of r->u and u->t, so this cannot form
                     * a cycle).  Keeping stale children made unpack explode. */
                    ex->via = u;
                    ex->child1 = e1;
                    ex->child2 = e2;
                }
            } else {
                uint32_t sid = ch_adj_add(a, uu, vv, len, t, fs, fc, fm, dx, dy,
                                          x->road_class, (uint8_t)(x->flags | y->flags), u,
                                          e1, e2);
                if (sid == VMR_INVALID) return -1;
                ch_pair_set(pairs, uu, vv, sid);
            }
            *(cur + uu) = (float)(*(a->out_deg + uu) + *(a->in_deg + uu));
            *(cur + vv) = (float)(*(a->out_deg + vv) + *(a->in_deg + vv));
            heap_push(pq, uu, *(cur + uu));
            heap_push(pq, vv, *(cur + vv));
        }
        ch_adj_kill(a, e1);
        *(cur + uu) = (float)(*(a->out_deg + uu) + *(a->in_deg + uu));
        heap_push(pq, uu, *(cur + uu));
    }
    for (uint32_t e2 = *(a->out_head + u); e2 != VMR_INVALID; e2 = *(a->out_next + e2)) {
        if (!(*(a->edges + e2)).alive) continue;
        uint32_t vv = (*(a->edges + e2)).to;
        ch_adj_kill(a, e2);
        *(cur + vv) = (float)(*(a->out_deg + vv) + *(a->in_deg + vv));
        heap_push(pq, vv, *(cur + vv));
    }
    *(contracted + u) = 1;
    *(level + u) = ++(*counter);
    return 0;
}

static int cmp_ch_packed(const void *pa, const void *pb) {
    const ChPacked *a = (const ChPacked *)pa;
    const ChPacked *b = (const ChPacked *)pb;
    if (a->from < b->from) return -1;
    if (a->from > b->from) return 1;
    if (a->to < b->to) return -1;
    if (a->to > b->to) return 1;
    return 0;
}

typedef struct { uint32_t from, to, old; } CchSortIdx;
static int cmp_cchsort(const void *pa, const void *pb) {
    const CchSortIdx *a = (const CchSortIdx *)pa;
    const CchSortIdx *b = (const CchSortIdx *)pb;
    if (a->from < b->from) return -1;
    if (a->from > b->from) return 1;
    if (a->to < b->to) return -1;
    if (a->to > b->to) return 1;
    return 0;
}

/* Vector-aware contraction hierarchy: shortcuts carry the same 8D features. */
static __attribute__((unused)) int vmr_build_ch_vector(VmrGraph *g) {
    uint32_t n = g->node_count;
    uint32_t m0 = g->edge_count;
    if (n == 0 || m0 == 0) {
        return -1;
    }
    ChAdj a;
    if (ch_adj_init(&a, n, m0) != 0) {
        ch_adj_free(&a);
        return -1;
    }
    for (uint32_t e = 0; e < m0; e++) {
        uint32_t from = *(g->edge_from + e);
        uint32_t to = *(g->edge_to + e);
        if (ch_adj_add(&a, from, to, *(g->edge_length_m + e), *(g->edge_time_s + e),
                       *(g->edge_f_speed + e), *(g->edge_f_class + e),
                       *(g->edge_f_smooth + e), *(g->edge_dx + e), *(g->edge_dy + e),
                       0, *(g->edge_flags + e), VMR_INVALID, VMR_INVALID, VMR_INVALID) == VMR_INVALID) {
            ch_adj_free(&a);
            return -1;
        }
    }

    uint8_t *contracted = (uint8_t *)xcalloc(n, sizeof(uint8_t));
    uint32_t *level = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint32_t *cur = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint8_t *delay = (uint8_t *)xcalloc(n, sizeof(uint8_t));
    ChPairSet pairs;
    MinHeap pq;
    if (!contracted || !level || !cur || !delay || ch_pair_init(&pairs, 1024) != 0 || heap_init(&pq, 1024) != 0) {
        xfree(contracted); xfree(level); xfree(cur); xfree(delay);
        ch_pair_free(&pairs);
        heap_free(&pq);
        ch_adj_free(&a);
        return -1;
    }

    for (uint32_t i = 0; i < n; i++) {
        *(cur + i) = (float)(*(a.out_deg + i) + *(a.in_deg + i));
        heap_push(&pq, i, *(cur + i));
    }
    LOG_INFO("CH: start (%u nodes, %u edges)", n, m0);

    uint32_t counter = 0;
    HeapEntry top;
    while (heap_pop(&pq, &top)) {
        uint32_t u = top.node;
        if (*(contracted + u) || top.prio != *(cur + u)) {
            continue;
        }
        uint32_t od = *(a.out_deg + u);
        uint32_t id = *(a.in_deg + u);
        if (od == 0 || id == 0) {
            if (ch_contract_node(&a, u, contracted, &pairs, level, &counter, cur, &pq) != 0) {
                goto fail;
            }
            continue;
        }
if ((uint64_t)od * (uint64_t)id > 32U && *(delay + u) < 8U) {
            (*(delay + u))++;
            *(cur + u) = (float)(od + id) * 1.5f;
            heap_push(&pq, u, *(cur + u));
            continue;
        }
        if (ch_contract_node(&a, u, contracted, &pairs, level, &counter, cur, &pq) != 0) {
            goto fail;
        }
            if ((counter & 65535U) == 0) {
                LOG_INFO("CH: %u contracted, %u edges", counter, (uint32_t)a.ecount);
            }
    }
    LOG_INFO("CH: queue empty, forcing remaining");
    for (uint32_t u = 0; u < n; u++) {
        if (!*(contracted + u)) {
            if (ch_contract_node(&a, u, contracted, &pairs, level, &counter, cur, &pq) != 0) {
                goto fail;
            }
            if ((counter & 65535U) == 0) {
                LOG_INFO("CH: forced %u contracted, %u edges", counter, (uint32_t)a.ecount);
            }
        }
    }
    LOG_INFO("CH: packing %u edges", (uint32_t)a.ecount);

    {
        size_t live = a.ecount; /* CH keeps every original edge and every shortcut */
        if (live == 0) {
            LOG_ERR("CH: no edges after contraction");
            goto fail;
        }
        LOG_INFO("CH: edges %u", (uint32_t)live);
        ChPacked *packed = (ChPacked *)xmalloc(live * sizeof(ChPacked));
        if (!packed) {
            LOG_ERR("CH: packed allocation failed");
            goto fail;
        }
        size_t w = 0;
        for (size_t e = 0; e < a.ecount; e++) {
            ChEdge *ce = a.edges + e;
            *(packed + w) = (ChPacked){
                ce->from, ce->to, ce->length_m, ce->time_s,
                ce->f_speed, ce->f_class, ce->f_smooth, ce->flags,
                ce->via, ce->child1, ce->child2
            };
            w++;
        }
        qsort(packed, live, sizeof(ChPacked), cmp_ch_packed);

        float *nlen = (float *)xmalloc(live * sizeof(float));
        float *ntime = (float *)xmalloc(live * sizeof(float));
        float *nfs = (float *)xmalloc(live * sizeof(float));
        float *nfc = (float *)xmalloc(live * sizeof(float));
        float *nfm = (float *)xmalloc(live * sizeof(float));
        float *ndx = (float *)xmalloc(live * sizeof(float));
        float *ndy = (float *)xmalloc(live * sizeof(float));
        float *nux = (float *)xmalloc(live * sizeof(float));
        float *nuy = (float *)xmalloc(live * sizeof(float));
        float *nfsa = (float *)xcalloc(live, sizeof(float));
        uint8_t *nfl = (uint8_t *)xmalloc(live);
        uint32_t *nfrom = (uint32_t *)xmalloc(live * sizeof(uint32_t));
        uint32_t *nto = (uint32_t *)xmalloc(live * sizeof(uint32_t));
        uint32_t *noff = (uint32_t *)xcalloc((size_t)n + 1U, sizeof(uint32_t));
        if (!nlen || !ntime || !nfs || !nfc || !nfm || !ndx || !ndy || !nux || !nuy ||
            !nfsa || !nfl || !nfrom || !nto || !noff) {
            LOG_ERR("CH: pack allocation failed (live=%u)", (uint32_t)live);
            xfree(packed); xfree(nlen); xfree(ntime); xfree(nfs); xfree(nfc); xfree(nfm);
            xfree(ndx); xfree(ndy); xfree(nux); xfree(nuy); xfree(nfsa); xfree(nfl);
            xfree(nfrom); xfree(nto); xfree(noff);
            goto fail;
        }
        for (size_t i = 0; i < live; i++) {
            ChPacked *p = packed + i;
            *(nfrom + i) = p->from;
            *(nto + i) = p->to;
            *(nlen + i) = p->length_m;
            *(ntime + i) = p->time_s;
            *(nfs + i) = p->f_speed;
            *(nfc + i) = p->f_class;
            *(nfm + i) = p->f_smooth;
            *(nfl + i) = p->flags;
            float dx = *(g->node_x + p->to) - *(g->node_x + p->from);
            float dy = *(g->node_y + p->to) - *(g->node_y + p->from);
            float mag = vmr_sqrt_fast(dx * dx + dy * dy);
            *(ndx + i) = dx;
            *(ndy + i) = dy;
            *(nux + i) = (mag > EPS_F) ? (dx / mag) : 0.0f;
            *(nuy + i) = (mag > EPS_F) ? (dy / mag) : 0.0f;
            float s = (*(nlen + i) > EPS_F) ? (mag / *(nlen + i)) : 1.0f;
            if (s > 1.0f) s = 1.0f;
            if (s < 0.05f) s = 0.05f;
            *(nfsa + i) = s;
            (*(noff + p->from + 1U))++;
        }
        for (uint32_t i = 0; i < n; i++) {
            *(noff + i + 1U) += *(noff + i);
        }

        if (!g->map) xfree(g->edge_from); xfree(g->edge_to);
        if (!g->map) xfree(g->edge_dx); xfree(g->edge_dy);
        if (!g->map) xfree(g->edge_ux); xfree(g->edge_uy);
        if (!g->map) xfree(g->edge_f_speed); xfree(g->edge_f_class);
        if (!g->map) xfree(g->edge_f_straight); xfree(g->edge_f_smooth);
        if (!g->map) xfree(g->edge_length_m); xfree(g->edge_time_s);
        if (!g->map) xfree(g->edge_flags); xfree(g->fwd_offsets);
        if (!g->map) xfree(g->rev_offsets); xfree(g->rev_src_node); xfree(g->rev_fwd_eid);

        g->edge_count = (uint32_t)live;
        g->edge_from = nfrom;
        g->edge_to = nto;
        g->edge_dx = ndx;
        g->edge_dy = ndy;
        g->edge_ux = nux;
        g->edge_uy = nuy;
        g->edge_f_speed = nfs;
        g->edge_f_class = nfc;
        g->edge_f_straight = nfsa;
        g->edge_f_smooth = nfm;
        g->edge_length_m = nlen;
        g->edge_time_s = ntime;
        g->edge_flags = nfl;
        g->fwd_offsets = noff;
        xfree(packed);
    }

    LOG_INFO("CH: packed %u edges, rebuilding reverse CSR", g->edge_count);
    if (vmr_build_reverse_csr_and_components(g) != 0) {
        LOG_ERR("CH: reverse CSR failed");
        goto fail;
    }
    LOG_INFO("CH: reverse CSR ok");
    g->ch_level = level;
    g->ch_built = 1;
    xfree(contracted);
    xfree(cur);
    xfree(delay);
    ch_pair_free(&pairs);
    heap_free(&pq);
    ch_adj_free(&a);
    LOG_INFO("CH built: %u nodes, %u edges (%u levels)",
             g->node_count, g->edge_count, counter);
    return 0;

fail:
    xfree(contracted); xfree(level); xfree(cur); xfree(delay);
    ch_pair_free(&pairs);
    heap_free(&pq);
    ch_adj_free(&a);
    return -1;
}

/* Scalar (travel-time) CH used only to build a coarse corridor.  The original
 * 8D vector graph is left untouched so the exact vector A* can refine inside
 * the corridor. */
static __attribute__((unused)) int vmr_build_ch(VmrGraph *g) {
    uint32_t n = g->node_count;
    uint32_t m0 = g->edge_count;
    if (n == 0 || m0 == 0) {
        return -1;
    }
    ChAdj a;
    if (ch_adj_init(&a, n, m0) != 0) {
        ch_adj_free(&a);
        return -1;
    }
    for (uint32_t e = 0; e < m0; e++) {
        uint32_t from = *(g->edge_from + e);
        uint32_t to = *(g->edge_to + e);
        if (ch_adj_add(&a, from, to, 0.0f, *(g->edge_time_s + e),
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0, 0, VMR_INVALID,
                       VMR_INVALID, VMR_INVALID) == VMR_INVALID) {
            ch_adj_free(&a);
            return -1;
        }
    }

    uint8_t *contracted = (uint8_t *)xcalloc(n, sizeof(uint8_t));
    uint32_t *level = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint32_t *cur = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint8_t *delay = (uint8_t *)xcalloc(n, sizeof(uint8_t));
    ChPairSet pairs;
    MinHeap pq;
    if (!contracted || !level || !cur || !delay ||
        ch_pair_init(&pairs, 1024) != 0 || heap_init(&pq, 1024) != 0) {
        xfree(contracted); xfree(level); xfree(cur); xfree(delay);
        ch_pair_free(&pairs); heap_free(&pq); ch_adj_free(&a);
        return -1;
    }

    for (uint32_t i = 0; i < n; i++) {
        *(cur + i) = (float)(*(a.out_deg + i) + *(a.in_deg + i));
        heap_push(&pq, i, *(cur + i));
    }

    uint32_t counter = 0;
    HeapEntry top;
    while (heap_pop(&pq, &top)) {
        heap_compact_min(&pq, n);
        uint32_t u = top.node;
        if (*(contracted + u) || top.prio != *(cur + u)) {
            continue;
        }
        uint32_t od = *(a.out_deg + u);
        uint32_t id = *(a.in_deg + u);
        if (od == 0 || id == 0) {
            if (ch_contract_node(&a, u, contracted, &pairs, level, &counter, cur, &pq) != 0) {
                goto fail;
            }
            continue;
        }
        if ((uint64_t)od * (uint64_t)id > 32U && *(delay + u) < 8U) {
            (*(delay + u))++;
            *(cur + u) = (float)(od + id) * 1.5f;
            heap_push(&pq, u, *(cur + u));
            continue;
        }
        if (ch_contract_node(&a, u, contracted, &pairs, level, &counter, cur, &pq) != 0) {
            goto fail;
        }
    }
    for (uint32_t u = 0; u < n; u++) {
        if (!*(contracted + u)) {
            if (ch_contract_node(&a, u, contracted, &pairs, level, &counter, cur, &pq) != 0) {
                goto fail;
            }
        }
    }

    {
        size_t total = a.ecount;
        if (total == 0) {
            goto fail;
        }
        /* Sort edges by (from,to) so the forward CSR groups each node's out-edges
         * contiguously, and REMAP every shortcut's child1/child2 through the same
         * permutation (otherwise children point at the wrong edge -> cycles). */
        CchSortIdx *sarr = (CchSortIdx *)xmalloc(total * sizeof(CchSortIdx));
        uint32_t *newpos = (uint32_t *)xmalloc(total * sizeof(uint32_t));
        if (!sarr || !newpos) {
            xfree(sarr); xfree(newpos);
            goto fail;
        }
        for (size_t e = 0; e < total; e++) {
            sarr[e].from = a.edges[e].from;
            sarr[e].to   = a.edges[e].to;
            sarr[e].old  = (uint32_t)e;
        }
        qsort(sarr, total, sizeof(CchSortIdx), cmp_cchsort);
        for (size_t j = 0; j < total; j++) {
            *(newpos + sarr[j].old) = (uint32_t)j;
        }

        uint32_t *cfrom = (uint32_t *)xmalloc(total * sizeof(uint32_t));
        uint32_t *cto = (uint32_t *)xmalloc(total * sizeof(uint32_t));
        float *ctime = (float *)xmalloc(total * sizeof(float));
        uint32_t *cvia = (uint32_t *)xmalloc(total * sizeof(uint32_t));
        uint32_t *cc1 = (uint32_t *)xmalloc(total * sizeof(uint32_t));
        uint32_t *cc2 = (uint32_t *)xmalloc(total * sizeof(uint32_t));
        uint32_t *coff = (uint32_t *)xcalloc((size_t)n + 1U, sizeof(uint32_t));
        if (!cfrom || !cto || !ctime || !cvia || !cc1 || !cc2 || !coff) {
            xfree(sarr); xfree(newpos);
            xfree(cfrom); xfree(cto); xfree(ctime); xfree(coff);
            xfree(cvia); xfree(cc1); xfree(cc2);
            goto fail;
        }
        for (size_t j = 0; j < total; j++) {
            ChEdge *ce = a.edges + sarr[j].old;
            *(cfrom + j) = ce->from;
            *(cto + j)   = ce->to;
            *(ctime + j) = ce->time_s;
            *(cvia + j)  = ce->via;
            *(cc1 + j)   = (ce->child1 == VMR_INVALID) ? VMR_INVALID : *(newpos + ce->child1);
            *(cc2 + j)   = (ce->child2 == VMR_INVALID) ? VMR_INVALID : *(newpos + ce->child2);
            (*(coff + ce->from + 1U))++;
        }
        xfree(sarr);
        xfree(newpos);
        for (uint32_t i = 0; i < n; i++) {
            *(coff + i + 1U) += *(coff + i);
        }

        uint32_t *roff = (uint32_t *)xcalloc((size_t)n + 1U, sizeof(uint32_t));
        uint32_t *rsrc = (uint32_t *)xmalloc(total * sizeof(uint32_t));
        uint32_t *reid = (uint32_t *)xmalloc(total * sizeof(uint32_t));
        uint32_t *rpos = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
        if (!roff || !rsrc || !reid || !rpos) {
            xfree(cfrom); xfree(cto); xfree(ctime); xfree(coff);
            xfree(cvia); xfree(cc1); xfree(cc2);
            xfree(roff); xfree(rsrc); xfree(reid); xfree(rpos);
            goto fail;
        }
        for (size_t i = 0; i < total; i++) {
            (*(roff + *(cto + i) + 1U))++;
        }
        for (uint32_t i = 0; i < n; i++) {
            *(roff + i + 1U) += *(roff + i);
        }
        memcpy(rpos, roff, (size_t)n * sizeof(uint32_t));
        for (size_t i = 0; i < total; i++) {
            uint32_t dst = *(cto + i);
            uint32_t at = (*(rpos + dst))++;
            *(rsrc + at) = *(cfrom + i);
            *(reid + at) = (uint32_t)i;
        }
        xfree(rpos);

        xfree(g->ch_from); xfree(g->ch_to); xfree(g->ch_time);
        xfree(g->ch_fwd_off); xfree(g->ch_rev_off);
        xfree(g->ch_rev_src); xfree(g->ch_rev_eid);
        xfree(g->ch_via); xfree(g->ch_child1); xfree(g->ch_child2);
        g->ch_edge_count = (uint32_t)total;
        g->ch_from = cfrom;
        g->ch_to = cto;
        g->ch_time = ctime;
        g->ch_fwd_off = coff;
        g->ch_rev_off = roff;
        g->ch_rev_src = rsrc;
        g->ch_rev_eid = reid;
        g->ch_via = cvia;
        g->ch_child1 = cc1;
        g->ch_child2 = cc2;
        g->ch_level = level;
        g->cch_built = 1;
    }

    xfree(contracted); xfree(cur); xfree(delay);
    ch_pair_free(&pairs); heap_free(&pq); ch_adj_free(&a);
    LOG_INFO("Scalar CH built: %u nodes, %u edges", g->node_count, g->ch_edge_count);
    return 0;

fail:
    xfree(contracted); xfree(level); xfree(cur); xfree(delay);
    ch_pair_free(&pairs); heap_free(&pq); ch_adj_free(&a);
    return -1;
}

/*
 * Performs full O(V + E) preprocessing on raw graph `rg` into compact SoA `g`:
 *   Step 1: Linked adjacency lists on raw edges for O(1) edge deletion/rewiring
 *   Step 2: Degree-2 chain contraction (one-way 1->m->1 and two-way 2<->m<->2)
 *   Step 3: Compaction of surviving nodes into dense indices
 *   Step 4: Dead-end 1-core cul-de-sac peeling
 *   Step 5: 8D SoA vector population, Reverse CSR, Components & Spatial Grid
 *
 * The raw graph storage is consumed and released before `g` is fully packed, so
 * preprocessing peak RSS stays close to a single graph copy instead of three.
 */
/* ------------------------------------------------------------------------- */
/* ALT landmarks (8 sources): exact landmark distances for goal-directed A*.  */
/* ------------------------------------------------------------------------- */

static uint32_t vmr_nearest_node_linear(const VmrGraph *g, double lat, double lon) {
    float kx = (float)(cos(lat * DEG_TO_RAD) * 111320.0);
    float ky = 110540.0f;
    float qx = (float)lon * kx;
    float qy = (float)lat * ky;
    float best = FLT_MAX;
    uint32_t bi = VMR_INVALID;
    for (uint32_t i = 0; i < g->node_count; i++) {
        float dx = (float)*(g->node_lon + i) * kx - qx;
        float dy = (float)*(g->node_lat + i) * ky - qy;
        float d2 = dx * dx + dy * dy;
        if (d2 < best) {
            best = d2;
            bi = i;
        }
    }
    return bi;
}

/* When building ALT landmarks we weight edges by time*attr_scale[0] so the
 * resulting heuristic is admissible (at scale 1.0) for the 8D cost
 * (time*attr*align*corridor >= time*attr). NULL => plain time. */
static const float *g_lm_attr = NULL;
static inline float alt_edge_w(const VmrGraph *g, uint32_t e) {
    float t = *(g->edge_time_s + e);
    return g_lm_attr ? t * g_lm_attr[e] : t;
}

static void alt_dijkstra(const VmrGraph *g, uint32_t src, float *dist, int reverse) {
    if (!dist || src >= g->node_count) {
        return;
    }
    uint32_t n = g->node_count;
    uint8_t *closed = (uint8_t *)xcalloc(n, 1);
    MinHeap pq;
    if (!closed || heap_init(&pq, 1024) != 0) {
        xfree(closed);
        return;
    }
    for (uint32_t i = 0; i < n; i++) {
        *(dist + i) = FLT_MAX;
    }
    *(dist + src) = 0.0f;
    heap_push(&pq, src, 0.0f);

    HeapEntry top;
    while (heap_pop(&pq, &top)) {
        uint32_t u = top.node;
        if (*(closed + u)) {
            continue;
        }
        *(closed + u) = 1;
        float du = *(dist + u);
        if (!reverse) {
            uint32_t beg = *(g->fwd_offsets + u);
            uint32_t end = *(g->fwd_offsets + u + 1U);
            for (uint32_t e = beg; e < end; e++) {
                uint32_t v = *(g->edge_to + e);
                if (*(closed + v)) {
                    continue;
                }
                float nd = du + alt_edge_w(g, e);
                if (nd < *(dist + v)) {
                    *(dist + v) = nd;
                    heap_push(&pq, v, nd);
                }
            }
        } else {
            uint32_t beg = *(g->rev_offsets + u);
            uint32_t end = *(g->rev_offsets + u + 1U);
            for (uint32_t r = beg; r < end; r++) {
                uint32_t v = *(g->rev_src_node + r);
                if (*(closed + v)) {
                    continue;
                }
                uint32_t e = *(g->rev_fwd_eid + r);
                float nd = du + alt_edge_w(g, e);
                if (nd < *(dist + v)) {
                    *(dist + v) = nd;
                    heap_push(&pq, v, nd);
                }
            }
        }
    }
    xfree(closed);
    heap_free(&pq);
}

static void vmr_build_landmarks(VmrGraph *g) {
    static const double lat[32] = {
        58.6440, 50.0700, 51.1300, 53.4800,
        55.9530, 51.5070, 52.4800, 57.4780,
        50.7200, 53.4000, 53.8000, 54.9700,
        55.8600, 57.1500, 51.4500, 52.6300,
        51.4816, 54.6032, 55.8617, 53.4808,
        52.4862, 53.8008, 53.4084, 54.9783,
        53.3811, 51.4545, 52.9548, 52.6369,
        50.8225, 57.1497, 57.4778, 51.6214
    };
    static const double lon[32] = {
        -3.0700, -5.7200, 1.3100, -2.2400,
        -3.1880, -0.1280, -1.9000, -4.2250,
        -3.5300, -2.9800, -1.5500, -1.6100,
        -4.2500, -2.0900, -2.5900, 1.3000,
        -3.1791, -5.9301, -4.2518, -2.2426,
        -1.8904, -1.5491, -2.9916, -1.6178,
        -1.4701, -2.5879, -1.1581, -1.1398,
        -0.1372, -2.0943, -4.2247, -3.9436
    };
    uint32_t nodes[32];
    uint32_t L = 0;
    for (uint32_t i = 0; i < 32; i++) {
        uint32_t u = vmr_nearest_node_linear(g, lat[i], lon[i]);
        if (u == VMR_INVALID) {
            continue;
        }
        int dup = 0;
        for (uint32_t j = 0; j < L; j++) {
            if (*(nodes + j) == u) {
                dup = 1;
                break;
            }
        }
        if (!dup) {
            *(nodes + L++) = u;
        }
    }
    if (L == 0) {
        return;
    }
    size_t n = g->node_count;
    g->landmark_node = (uint32_t *)xmalloc((size_t)L * sizeof(uint32_t));
    g->alt_nodes = (AltNode32 *)xaligned_alloc(64, n * sizeof(AltNode32));
    if (!g->landmark_node || !g->alt_nodes) {
        if (!g->map) xfree(g->landmark_node);
        xfree_aligned(g->alt_nodes);
        g->landmark_node = NULL;
        g->alt_nodes = NULL;
        g->landmark_count = 0;
        LOG_WARN("ALT landmark allocation failed; falling back to Euclidean heuristic");
        return;
    }
    memcpy(g->landmark_node, nodes, (size_t)L * sizeof(uint32_t));
    for (size_t v = 0; v < n; v++) {
        for (int k = 0; k < 32; k++) {
            (*(g->alt_nodes + v)).fwd[k] = 1.0e30f;
            (*(g->alt_nodes + v)).bwd[k] = -1.0e30f;
        }
    }
    g->landmark_count = L;
    /* Time-based landmarks (admissible for all modes). */
    g_lm_attr = NULL;
    clock_t t0 = clock();
#pragma omp parallel for schedule(dynamic)
    for (int li = 0; li < (int)L; li++) {
        uint32_t u = *(g->landmark_node + li);
        float *tmp = (float *)xmalloc(n * sizeof(float));
        if (!tmp) {
            continue;
        }
        alt_dijkstra(g, u, tmp, 0);
        for (size_t v = 0; v < n; v++) {
            float d = *(tmp + v);
            (*(g->alt_nodes + v)).fwd[li] = (d >= 1.0e29f) ? 1.0e30f : d;
        }
        alt_dijkstra(g, u, tmp, 1);
        for (size_t v = 0; v < n; v++) {
            float d = *(tmp + v);
            (*(g->alt_nodes + v)).bwd[li] = (d >= 1.0e29f) ? -1.0e30f : d;
        }
        xfree(tmp);
    }
    clock_t t1 = clock();
    LOG_INFO("ALT landmarks: %u built in %.1f ms", L,
             1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC);
    g_lm_attr = NULL;
}

static void vmr_build_edge_blocks(VmrGraph *g) {
    uint32_t n = g->node_count;
    size_t rnb = 0;
    for (uint32_t w = 0; w < n; w++) {
        uint32_t deg = *(g->rev_offsets + w + 1U) - *(g->rev_offsets + w);
        rnb += (deg + 3U) / 4U;
    }
    EdgeBlock4 *rblocks = (EdgeBlock4 *)xaligned_alloc(32, (rnb ? rnb : 1U) * sizeof(EdgeBlock4));
    uint32_t *rboff = (uint32_t *)xcalloc((size_t)n + 1U, sizeof(uint32_t));
    if (rblocks && rboff) {
        memset(rblocks, 0, (rnb ? rnb : 1U) * sizeof(EdgeBlock4));
        size_t rb = 0;
        for (uint32_t w = 0; w < n; w++) {
            *(rboff + w) = (uint32_t)rb;
            uint32_t beg = *(g->rev_offsets + w);
            uint32_t end = *(g->rev_offsets + w + 1U);
            for (uint32_t ri = beg; ri < end; ri += 4U) {
                EdgeBlock4 *blk = rblocks + rb;
                for (uint32_t l = 0; l < 4; l++) {
                    if (ri + l < end) {
                        uint32_t eid = *(g->rev_fwd_eid + ri + l);
                        uint32_t to = *(g->rev_src_node + ri + l);
                        blk->to[l] = to;
                        blk->ux[l] = -(*(g->edge_ux + eid));
                        blk->uy[l] = -(*(g->edge_uy + eid));
                        blk->dx[l] = -(*(g->edge_dx + eid));
                        blk->dy[l] = -(*(g->edge_dy + eid));
                        blk->vx[l] = *(g->node_x + to);
                        blk->vy[l] = *(g->node_y + to);
                        blk->time[l] = *(g->edge_time_s + eid);
                        for (uint32_t p = 0; p < 5; p++) {
                            blk->attr[p][l] = g->edge_attr_scale[p]
                                                  ? *(g->edge_attr_scale[p] + eid) : 1.0f;
                        }
                        blk->flags[l] = *(g->edge_flags + eid);
                    } else {
                        blk->to[l] = VMR_INVALID;
                        blk->flags[l] = EDGE_DEAD_SOURCE;
                    }
                }
                rb++;
            }
        }
        *(rboff + n) = (uint32_t)rb;
        xfree_aligned(g->rev_blocks);
        xfree(g->rev_block_off);
        g->rev_blocks = rblocks;
        g->rev_block_off = rboff;
    } else {
        xfree_aligned(rblocks);
        xfree(rboff);
    }

    LOG_INFO("Edge blocks: %u rev (%u edges)", (uint32_t)rnb, g->edge_count);
}

static void vmr_build_two_hop_rev(VmrGraph *g) {
    size_t m = g->edge_count;
    if (m == 0) {
        return;
    }
    uint32_t *off = (uint32_t *)xmalloc((m + 1U) * sizeof(uint32_t));
    if (!off) {
        return;
    }
    /* Pass 1: count continuations per reverse edge (parallel). */
    #pragma omp parallel for schedule(static)
    for (uint32_t ri = 0; ri < m; ri++) {
        uint32_t eid = *(g->rev_fwd_eid + ri);
        uint32_t w = *(g->edge_from + eid);
        uint32_t u = *(g->rev_src_node + ri);
        uint32_t beg = *(g->rev_offsets + u);
        uint32_t end = *(g->rev_offsets + u + 1U);
        uint32_t cnt = 0;
        if (end - beg <= 8U) {
            const uint32_t *rs = g->rev_src_node;
            for (uint32_t r2 = beg; r2 < end; r2++) {
                if (*(rs + r2) != w) {
                    cnt++;
                }
            }
        }
        *(off + ri + 1U) = cnt;
    }
    uint32_t total = 0;
    for (uint32_t ri = 0; ri < m; ri++) {
        uint32_t cnt = *(off + ri + 1U);
        *(off + ri) = total;
        total += cnt;
    }
    *(off + m) = total;
    float *sx = (float *)xmalloc((size_t)(total ? total : 1U) * sizeof(float));
    float *sy = (float *)xmalloc((size_t)(total ? total : 1U) * sizeof(float));
    float *iv = (float *)xmalloc((size_t)(total ? total : 1U) * sizeof(float));
    if (!sx || !sy || !iv) {
        xfree(off); xfree(sx); xfree(sy); xfree(iv);
        return;
    }
    /* Pass 2: fill (parallel). */
    #pragma omp parallel for schedule(static)
    for (uint32_t ri = 0; ri < m; ri++) {
        uint32_t c = *(off + ri);
        uint32_t eid = *(g->rev_fwd_eid + ri);
        uint32_t w = *(g->edge_from + eid);
        uint32_t u = *(g->rev_src_node + ri);
        uint32_t beg = *(g->rev_offsets + u);
        uint32_t end = *(g->rev_offsets + u + 1U);
        if (end - beg > 8U) {
            continue;
        }
        float d1x = -(*(g->edge_dx + eid));
        float d1y = -(*(g->edge_dy + eid));
        const uint32_t *rs = g->rev_src_node;
        const uint32_t *re = g->rev_fwd_eid;
        const float *dxp = g->edge_dx, *dyp = g->edge_dy;
        for (uint32_t r2 = beg; r2 < end; r2++) {
            uint32_t prev = *(rs + r2);
            if (prev == w) continue;
            uint32_t eid2 = *(re + r2);
            float ax = d1x - *(dxp + eid2);
            float ay = d1y - *(dyp + eid2);
            *(sx + c) = ax;
            *(sy + c) = ay;
            *(iv + c) = 1.0f / (sqrtf(ax * ax + ay * ay) + EPS_F);
            c++;
        }
    }
    xfree(g->rth_off);
    xfree(g->rth_sx);
    xfree(g->rth_sy);
    xfree(g->rth_inv);
    g->rth_off = off;
    g->rth_sx = sx;
    g->rth_sy = sy;
    g->rth_inv = iv;
    g->rth_count = total;
    LOG_INFO("Two-hop reverse continuations: %u entries", total);
}

/*
 * Precompute the 2-hop continuation vectors for every edge. The vector sum
 * (dx1+dx2, dy1+dy2) and its inverse norm are query-independent, so the hot
 * loop only performs a dot product instead of an adjacency scan + sqrt.
 */
static void vmr_build_two_hop(VmrGraph *g) {
    size_t m = g->edge_count;
    if (m == 0) {
        return;
    }
    uint32_t *off = (uint32_t *)xmalloc((m + 1U) * sizeof(uint32_t));
    if (!off) {
        return;
    }
    /* Pass 1: count continuations per edge (parallel). */
    #pragma omp parallel for schedule(static)
    for (uint32_t e = 0; e < m; e++) {
        uint32_t v = *(g->edge_to + e);
        uint32_t beg = *(g->fwd_offsets + v);
        uint32_t end = *(g->fwd_offsets + v + 1U);
        uint32_t cnt = 0;
        if (end - beg <= 8U) {
            uint32_t u = *(g->edge_from + e);
            const uint32_t *to_p = g->edge_to;
            for (uint32_t e2 = beg; e2 < end; e2++) {
                if (*(to_p + e2) != u) {
                    cnt++;
                }
            }
        }
        *(off + e + 1U) = cnt;
    }
    /* Prefix sum -> start offset per edge. */
    uint32_t total = 0;
    for (uint32_t e = 0; e < m; e++) {
        uint32_t cnt = *(off + e + 1U);
        *(off + e) = total;
        total += cnt;
    }
    *(off + m) = total;
    float *sx = (float *)xmalloc((size_t)(total ? total : 1U) * sizeof(float));
    float *sy = (float *)xmalloc((size_t)(total ? total : 1U) * sizeof(float));
    float *iv = (float *)xmalloc((size_t)(total ? total : 1U) * sizeof(float));
    if (!sx || !sy || !iv) {
        xfree(off); xfree(sx); xfree(sy); xfree(iv);
        return;
    }
    /* Pass 2: fill (parallel). */
    #pragma omp parallel for schedule(static)
    for (uint32_t e = 0; e < m; e++) {
        uint32_t c = *(off + e);
        uint32_t u = *(g->edge_from + e);
        uint32_t v = *(g->edge_to + e);
        uint32_t beg = *(g->fwd_offsets + v);
        uint32_t end = *(g->fwd_offsets + v + 1U);
        if (end - beg > 8U) {
            continue;
        }
        float d1x = *(g->edge_dx + e), d1y = *(g->edge_dy + e);
        const uint32_t *to_p = g->edge_to;
        const float *dx_p = g->edge_dx, *dy_p = g->edge_dy;
        for (uint32_t e2 = beg; e2 < end; e2++) {
            uint32_t w = *(to_p + e2);
            if (w == u) continue;
            float ax = d1x + *(dx_p + e2);
            float ay = d1y + *(dy_p + e2);
            *(sx + c) = ax;
            *(sy + c) = ay;
            *(iv + c) = 1.0f / (sqrtf(ax * ax + ay * ay) + EPS_F);
            c++;
        }
    }
    xfree(g->th_off);
    xfree(g->th_sx);
    xfree(g->th_sy);
    xfree(g->th_inv);
    g->th_off = off;
    g->th_sx = sx;
    g->th_sy = sy;
    g->th_inv = iv;
    g->th_count = total;
    LOG_INFO("Two-hop continuations: %u entries", total);
}

typedef struct {
    uint32_t id;
    uint32_t d;
} VmrHNode;

static int vmr_hnode_cmp(const void *a, const void *b) {
    const VmrHNode *x = (const VmrHNode *)a;
    const VmrHNode *y = (const VmrHNode *)b;
    if (x->d != y->d) {
        return (x->d < y->d) ? -1 : 1;
    }
    return (x->id > y->id) - (x->id < y->id);
}

static uint32_t vmr_hilbert_xy2d(uint32_t x, uint32_t y, uint32_t order) {
    uint32_t n = 1u << order;
    uint32_t d = 0;
    for (uint32_t s = n >> 1; s; s >>= 1) {
        uint32_t rx = (x & s) ? 1u : 0u;
        uint32_t ry = (y & s) ? 1u : 0u;
        d += s * s * ((3u * rx) ^ ry);
        if (ry == 0) {
            if (rx) {
                x = n - 1u - x;
                y = n - 1u - y;
            }
            uint32_t t = x;
            x = y;
            y = t;
        }
    }
    return d;
}

/*
 * Locality reorder: assign new node ids along a Hilbert curve over the
 * projected (x, y) plane, and permute every per-node/per-edge array to match.
 * The search frontier then touches nearby cache lines, which is the dominant
 * cost at 84k expansions on a 2.7M-node graph.
 */
static void vmr_reorder_graph(VmrGraph *g) {
    uint32_t n = g->node_count;
    size_t m = g->edge_count;
    if (n < 2 || m == 0 || !g->fwd_offsets || !g->node_x || !g->node_y) {
        return;
    }
    float minx = *(g->node_x), maxx = minx;
    float miny = *(g->node_y), maxy = miny;
    for (uint32_t i = 1; i < n; i++) {
        float x = *(g->node_x + i), y = *(g->node_y + i);
        if (x < minx) minx = x;
        if (x > maxx) maxx = x;
        if (y < miny) miny = y;
        if (y > maxy) maxy = y;
    }
    float scx = (maxx > minx) ? 65535.0f / (maxx - minx) : 0.0f;
    float scy = (maxy > miny) ? 65535.0f / (maxy - miny) : 0.0f;

    VmrHNode *hn = (VmrHNode *)xmalloc((size_t)n * sizeof(VmrHNode));
    uint32_t *perm = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint32_t *iperm = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint32_t *nfill = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint32_t *rfill = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    if (!hn || !perm || !iperm || !nfill || !rfill) {
        xfree(hn); xfree(perm); xfree(iperm); xfree(nfill); xfree(rfill);
        return;
    }
    for (uint32_t i = 0; i < n; i++) {
        hn[i].id = i;
        uint32_t hx = (uint32_t)((*(g->node_x + i) - minx) * scx);
        uint32_t hy = (uint32_t)((*(g->node_y + i) - miny) * scy);
        hn[i].d = vmr_hilbert_xy2d(hx, hy, 16);
    }
    qsort(hn, n, sizeof(VmrHNode), vmr_hnode_cmp);
    for (uint32_t i = 0; i < n; i++) {
        uint32_t old = hn[i].id;
        perm[old] = i;
        iperm[i] = old;
    }
    xfree(hn);

    int64_t *n_oid = (int64_t *)xmalloc((size_t)n * sizeof(int64_t));
    double *n_lat = (double *)xmalloc((size_t)n * sizeof(double));
    double *n_lon = (double *)xmalloc((size_t)n * sizeof(double));
    float *n_x = (float *)xmalloc((size_t)n * sizeof(float));
    float *n_y = (float *)xmalloc((size_t)n * sizeof(float));
    uint32_t *n_comp = (uint32_t *)xmalloc((size_t)n * sizeof(uint32_t));
    uint8_t *n_flags = (uint8_t *)xmalloc((size_t)n);
    if (!n_oid || !n_lat || !n_lon || !n_x || !n_y || !n_comp || !n_flags) {
        xfree(n_oid); xfree(n_lat); xfree(n_lon); xfree(n_x); xfree(n_y);
        xfree(n_comp); xfree(n_flags);
        xfree(perm); xfree(iperm); xfree(nfill); xfree(rfill);
        return;
    }
    for (uint32_t old = 0; old < n; old++) {
        uint32_t i = perm[old];
        *(n_oid + i) = *(g->node_osm_id + old);
        *(n_lat + i) = *(g->node_lat + old);
        *(n_lon + i) = *(g->node_lon + old);
        *(n_x + i) = *(g->node_x + old);
        *(n_y + i) = *(g->node_y + old);
        *(n_comp + i) = *(g->node_comp_id + old);
        *(n_flags + i) = *(g->node_flags + old);
    }

    uint32_t *noff = (uint32_t *)xcalloc((size_t)n + 1U, sizeof(uint32_t));
    if (!noff) {
        xfree(n_oid); xfree(n_lat); xfree(n_lon); xfree(n_x); xfree(n_y);
        xfree(n_comp); xfree(n_flags);
        xfree(perm); xfree(iperm); xfree(nfill); xfree(rfill);
        return;
    }
    for (uint32_t old = 0; old < n; old++) {
        uint32_t deg = *(g->fwd_offsets + old + 1U) - *(g->fwd_offsets + old);
        *(noff + perm[old] + 1U) = deg;
    }
    for (uint32_t i = 0; i < n; i++) {
        *(noff + i + 1U) += *(noff + i);
    }
    memcpy(nfill, noff, (size_t)n * sizeof(uint32_t));

    uint32_t *n_from = (uint32_t *)xmalloc(m * sizeof(uint32_t));
    uint32_t *n_to = (uint32_t *)xmalloc(m * sizeof(uint32_t));
    float *n_dx = (float *)xmalloc(m * sizeof(float));
    float *n_dy = (float *)xmalloc(m * sizeof(float));
    float *n_ux = (float *)xmalloc(m * sizeof(float));
    float *n_uy = (float *)xmalloc(m * sizeof(float));
    float *n_fs = (float *)xmalloc(m * sizeof(float));
    float *n_fc = (float *)xmalloc(m * sizeof(float));
    float *n_fst = (float *)xmalloc(m * sizeof(float));
    float *n_fsm = (float *)xmalloc(m * sizeof(float));
    float *n_len = (float *)xmalloc(m * sizeof(float));
    float *n_time = (float *)xmalloc(m * sizeof(float));
    uint8_t *n_eflags = (uint8_t *)xmalloc(m);
    float *n_attr[5];
    for (uint32_t p = 0; p < 5; p++) {
        n_attr[p] = (float *)xmalloc(m * sizeof(float));
    }
    if (!n_from || !n_to || !n_dx || !n_dy || !n_ux || !n_uy || !n_fs || !n_fc ||
        !n_fst || !n_fsm || !n_len || !n_time || !n_eflags) {
        xfree(n_from); xfree(n_to); xfree(n_dx); xfree(n_dy); xfree(n_ux); xfree(n_uy);
        xfree(n_fs); xfree(n_fc); xfree(n_fst); xfree(n_fsm); xfree(n_len);
        xfree(n_time); xfree(n_eflags);
        xfree(n_oid); xfree(n_lat); xfree(n_lon); xfree(n_x); xfree(n_y);
        xfree(n_comp); xfree(n_flags); xfree(noff);
        xfree(perm); xfree(iperm); xfree(nfill); xfree(rfill);
        return;
    }
    for (uint32_t p = 0; p < 5; p++) {
        if (n_attr[p]) {
            continue;
        }
        xfree(n_from); xfree(n_to); xfree(n_dx); xfree(n_dy); xfree(n_ux); xfree(n_uy);
        xfree(n_fs); xfree(n_fc); xfree(n_fst); xfree(n_fsm); xfree(n_len);
        xfree(n_time); xfree(n_eflags);
        xfree(n_oid); xfree(n_lat); xfree(n_lon); xfree(n_x); xfree(n_y);
        xfree(n_comp); xfree(n_flags); xfree(noff);
        xfree(perm); xfree(iperm); xfree(nfill); xfree(rfill);
        return;
    }
    for (uint32_t old = 0; old < n; old++) {
        uint32_t u2 = perm[old];
        uint32_t beg = *(g->fwd_offsets + old);
        uint32_t end = *(g->fwd_offsets + old + 1U);
        for (uint32_t e = beg; e < end; e++) {
            uint32_t pos = *(nfill + u2);
            *(nfill + u2) = pos + 1U;
            *(n_from + pos) = u2;
            *(n_to + pos) = perm[*(g->edge_to + e)];
            *(n_dx + pos) = *(g->edge_dx + e);
            *(n_dy + pos) = *(g->edge_dy + e);
            *(n_ux + pos) = *(g->edge_ux + e);
            *(n_uy + pos) = *(g->edge_uy + e);
            *(n_fs + pos) = *(g->edge_f_speed + e);
            *(n_fc + pos) = *(g->edge_f_class + e);
            *(n_fst + pos) = *(g->edge_f_straight + e);
            *(n_fsm + pos) = *(g->edge_f_smooth + e);
            *(n_len + pos) = *(g->edge_length_m + e);
            *(n_time + pos) = *(g->edge_time_s + e);
            *(n_eflags + pos) = *(g->edge_flags + e);
            for (uint32_t p = 0; p < 5; p++) {
                *(n_attr[p] + pos) = g->edge_attr_scale[p]
                                         ? *(g->edge_attr_scale[p] + e) : 1.0f;
            }
        }
    }

    uint32_t *roff = (uint32_t *)xcalloc((size_t)n + 1U, sizeof(uint32_t));
    uint32_t *rsrc = (uint32_t *)xmalloc(m * sizeof(uint32_t));
    uint32_t *reid = (uint32_t *)xmalloc(m * sizeof(uint32_t));
    if (roff && rsrc && reid) {
        for (size_t e = 0; e < m; e++) {
            *(roff + *(n_to + e) + 1U) += 1U;
        }
        for (uint32_t i = 0; i < n; i++) {
            *(roff + i + 1U) += *(roff + i);
        }
        memcpy(rfill, roff, (size_t)n * sizeof(uint32_t));
        for (size_t e = 0; e < m; e++) {
            uint32_t v = *(n_to + e);
            uint32_t r = *(rfill + v);
            *(rfill + v) = r + 1U;
            *(rsrc + r) = *(n_from + e);
            *(reid + r) = (uint32_t)e;
        }
    } else {
        xfree(roff);
        roff = NULL;
        xfree(rsrc);
        xfree(reid);
    }

    uint32_t *ntoff = NULL;
    uint32_t *ntfrom = NULL;
    uint32_t *ntto = NULL;
    uint8_t *ntmode = NULL;
    if (g->turn_count > 0 && g->turn_offsets && g->turn_from && g->turn_to && g->turn_mode) {
        ntoff = (uint32_t *)xmalloc((size_t)(n + 1U) * sizeof(uint32_t));
        ntfrom = (uint32_t *)xmalloc(g->turn_count * sizeof(uint32_t));
        ntto = (uint32_t *)xmalloc(g->turn_count * sizeof(uint32_t));
        ntmode = (uint8_t *)xmalloc(g->turn_count);
        if (ntoff && ntfrom && ntto && ntmode) {
            uint32_t tc = 0;
            for (uint32_t v = 0; v < n; v++) {
                uint32_t old = iperm[v];
                *(ntoff + v) = tc;
                uint32_t beg = *(g->turn_offsets + old);
                uint32_t end = *(g->turn_offsets + old + 1U);
                for (uint32_t i = beg; i < end; i++) {
                    *(ntfrom + tc) = perm[*(g->turn_from + i)];
                    *(ntto + tc) = perm[*(g->turn_to + i)];
                    *(ntmode + tc) = *(g->turn_mode + i);
                    tc++;
                }
            }
            *(ntoff + n) = tc;
        } else {
            xfree(ntoff);
            ntoff = NULL;
            xfree(ntfrom);
            ntfrom = NULL;
            xfree(ntto);
            ntto = NULL;
            xfree(ntmode);
            ntmode = NULL;
        }
    }

    if (!g->map) xfree(g->node_osm_id);
    if (!g->map) xfree(g->node_lat);
    if (!g->map) xfree(g->node_lon);
    if (!g->map) xfree(g->node_x);
    if (!g->map) xfree(g->node_y);
    if (!g->map) xfree(g->node_comp_id);
    if (!g->map) xfree(g->node_flags);
    if (!g->map) xfree(g->fwd_offsets);
    if (!g->map) xfree(g->edge_from);
    if (!g->map) xfree(g->edge_to);
    if (!g->map) xfree(g->edge_dx);
    if (!g->map) xfree(g->edge_dy);
    if (!g->map) xfree(g->edge_ux);
    if (!g->map) xfree(g->edge_uy);
    if (!g->map) xfree(g->edge_f_speed);
    if (!g->map) xfree(g->edge_f_class);
    if (!g->map) xfree(g->edge_f_straight);
    if (!g->map) xfree(g->edge_f_smooth);
    if (!g->map) xfree(g->edge_length_m);
    if (!g->map) xfree(g->edge_time_s);
    if (!g->map) xfree(g->edge_flags);
    for (uint32_t p = 0; p < 5; p++) {
        if (!g->map) xfree(g->edge_attr_scale[p]);
    }
    if (!g->map) xfree(g->rev_offsets);
    if (!g->map) xfree(g->rev_src_node);
    if (!g->map) xfree(g->rev_fwd_eid);
    xfree(g->turn_offsets);
    xfree(g->turn_from);
    xfree(g->turn_to);
    xfree(g->turn_mode);

    g->node_osm_id = n_oid;
    g->node_lat = n_lat;
    g->node_lon = n_lon;
    g->node_x = n_x;
    g->node_y = n_y;
    g->node_comp_id = n_comp;
    g->node_flags = n_flags;
    g->fwd_offsets = noff;
    g->edge_from = n_from;
    g->edge_to = n_to;
    g->edge_dx = n_dx;
    g->edge_dy = n_dy;
    g->edge_ux = n_ux;
    g->edge_uy = n_uy;
    g->edge_f_speed = n_fs;
    g->edge_f_class = n_fc;
    g->edge_f_straight = n_fst;
    g->edge_f_smooth = n_fsm;
    g->edge_length_m = n_len;
    g->edge_time_s = n_time;
    g->edge_flags = n_eflags;
    for (uint32_t p = 0; p < 5; p++) {
        g->edge_attr_scale[p] = n_attr[p];
    }
    g->rev_offsets = roff;
    g->rev_src_node = rsrc;
    g->rev_fwd_eid = reid;
    g->turn_offsets = ntoff;
    g->turn_from = ntfrom;
    g->turn_to = ntto;
    g->turn_mode = ntmode;

    xfree(perm);
    xfree(iperm);
    xfree(nfill);
    xfree(rfill);
    LOG_INFO("Hilbert reorder: %u nodes, %u edges", n, (uint32_t)m);
}

static int vmr_blocks_enabled(void) {
    static int v = -1;
    if (v < 0) {
        v = (getenv("VMR_BLOCK") != NULL) ? 1 : 0;
    }
    return v;
}

static int vmr_rev_enabled(void) {
    static int v = -1;
    if (v < 0) {
        const char *s = getenv("VMR_BLOCK");
        v = (s && (strcmp(s, "rev") == 0 || strcmp(s, "2") == 0 ||
                   strcmp(s, "full") == 0)) ? 1 : 0;
    }
    return v;
}

static int vmr_preprocess_graph(VmrRawGraph *rg, VmrGraph *g) {
    vmr_graph_free(g);
    vmr_graph_init(g);

    uint32_t rn = rg->node_count;
    uint32_t rm = rg->edge_count;
    if (rn == 0 || rm == 0) {
        LOG_ERR("Cannot preprocess empty graph");
        return -1;
    }

    clock_t t0 = clock();
    g->raw_node_count = rn;
    g->raw_edge_count = rm;

    /* Contract in place on the raw edge pool: no second 32 B/edge copy. */
    RawEdge *we       = rg->edges;
    uint32_t *out_head = (uint32_t *)xmalloc((size_t)rn * sizeof(uint32_t));
    uint32_t *in_head  = (uint32_t *)xmalloc((size_t)rn * sizeof(uint32_t));
    uint32_t *out_next = (uint32_t *)xmalloc((size_t)rm * sizeof(uint32_t));
    uint32_t *in_next  = (uint32_t *)xmalloc((size_t)rm * sizeof(uint32_t));
    uint32_t *out_deg  = (uint32_t *)xcalloc(rn, sizeof(uint32_t));
    uint32_t *in_deg   = (uint32_t *)xcalloc(rn, sizeof(uint32_t));
    uint8_t  *removed  = (uint8_t *)xcalloc(rn, sizeof(uint8_t));

    if (!we || !out_head || !in_head || !out_next || !in_next ||
        !out_deg || !in_deg || !removed) {
        xfree(out_head); xfree(in_head); xfree(out_next);
        xfree(in_next); xfree(out_deg); xfree(in_deg); xfree(removed);
        return -1;
    }

    for (uint32_t i = 0; i < rn; i++) {
        *(out_head + i) = VMR_INVALID;
        *(in_head + i)  = VMR_INVALID;
    }

    for (uint32_t e = 0; e < rm; e++) {
        RawEdge *w = we + e;
        w->alive = 1;

        *(out_next + e) = *(out_head + w->from);
        *(out_head + w->from) = e;
        (*(out_deg + w->from))++;

        *(in_next + e) = *(in_head + w->to);
        *(in_head + w->to) = e;
        (*(in_deg + w->to))++;
    }

    /* Step 2: Degree-2 Chain Contraction (iterative pass until fixed point) */
    uint32_t contracted_cnt = 0;
    uint8_t *no_contract = NULL;
    if (rg->turn_count > 0) {
        no_contract = (uint8_t *)xcalloc(rn, sizeof(uint8_t));
        if (no_contract) {
            for (uint32_t i = 0; i < rg->turn_count; i++) {
                *(no_contract + *(rg->turn_via + i)) = 1;
            }
        }
    }
    uint32_t *cq = (uint32_t *)xmalloc((size_t)rn * sizeof(uint32_t));
    uint8_t  *in_cq = (uint8_t *)xcalloc(rn, sizeof(uint8_t));
    if (!cq || !in_cq) {
        xfree(cq); xfree(in_cq);
        xfree(out_head); xfree(in_head); xfree(out_next);
        xfree(in_next); xfree(out_deg); xfree(in_deg); xfree(removed);
        return -1;
    }

    uint32_t cq_head = 0, cq_tail = 0;
    for (uint32_t i = 0; i < rn; i++) {
        uint32_t od = *(out_deg + i);
        uint32_t id = *(in_deg + i);
        if ((od == 1 && id == 1) || (od == 2 && id == 2)) {
            if (no_contract && *(no_contract + i)) continue;
            *(cq + cq_tail++) = i;
            *(in_cq + i) = 1;
        }
    }

    while (cq_head < cq_tail) {
        uint32_t m = *(cq + cq_head++);
        *(in_cq + m) = 0;
        if (*(removed + m)) {
            continue;
        }
        uint32_t od = *(out_deg + m);
        uint32_t id = *(in_deg + m);

        /* Collect live incoming and outgoing edge indices for node m */
        uint32_t in_e[2]  = { VMR_INVALID, VMR_INVALID };
        uint32_t out_e[2] = { VMR_INVALID, VMR_INVALID };
        uint32_t ic = 0, oc = 0;

        for (uint32_t e = *(in_head + m); e != VMR_INVALID && ic < 3; e = *(in_next + e)) {
            if ((we + e)->alive) {
                if (ic < 2) in_e[ic] = e;
                ic++;
            }
        }
        for (uint32_t e = *(out_head + m); e != VMR_INVALID && oc < 3; e = *(out_next + e)) {
            if ((we + e)->alive) {
                if (oc < 2) out_e[oc] = e;
                oc++;
            }
        }

        if (ic == 1 && oc == 1) {
            /* One-way chain: u -> m -> v */
            uint32_t e1 = in_e[0], e2 = out_e[0];
            RawEdge *w1 = we + e1;
            RawEdge *w2 = we + e2;
            uint32_t u = w1->from;
            uint32_t v = w2->to;
            if (u == v || u == m || v == m) continue;
            if (w1->road_class != w2->road_class ||
                (w1->flags & EDGE_TOLL) != (w2->flags & EDGE_TOLL)) {
                continue;
            }
            /* Merge w1 + w2 into w1 (now u -> v), kill w2 */
            float tot_len = w1->length_m + w2->length_m;
            w1->f_speed   = (w1->f_speed * w1->length_m + w2->f_speed * w2->length_m) / tot_len;
            w1->f_smooth  = (w1->f_smooth * w1->length_m + w2->f_smooth * w2->length_m) / tot_len;
            w1->length_m  = tot_len;
            w1->travel_time_s += w2->travel_time_s;
            w1->to        = v;
            w1->flags    |= EDGE_CONTRACTED | (w2->flags & EDGE_BRIDGE);

            /* Rewire w1 from m's in-list into v's in-list */
            *(in_next + e1) = *(in_head + v);
            *(in_head + v)  = e1;
            /* Note: in_deg[v] had +1 from e2 and loses e2, gains e1 -> unchanged! */

            w2->alive = 0;
            *(out_deg + m) = 0;
            *(in_deg + m)  = 0;
            *(removed + m) = 1;
            contracted_cnt++;
        } else if (ic == 2 && oc == 2 && od == 2 && id == 2) {
            /* Two-way chain: u <-> m <-> v */
            RawEdge *wi0 = we + in_e[0];
            RawEdge *wi1 = we + in_e[1];
            RawEdge *wo0 = we + out_e[0];
            RawEdge *wo1 = we + out_e[1];

            uint32_t u = wi0->from;
            uint32_t v = wi1->from;
            if (u == v || u == m || v == m) continue;

            /* Match outgoing edges so e_mv goes to v and e_mu goes to u */
            uint32_t e_um = in_e[0], e_vm = in_e[1];
            uint32_t e_mv = VMR_INVALID, e_mu = VMR_INVALID;
            if (wo0->to == v && wo1->to == u) {
                e_mv = out_e[0];
                e_mu = out_e[1];
            } else if (wo0->to == u && wo1->to == v) {
                e_mu = out_e[0];
                e_mv = out_e[1];
            } else {
                continue;
            }

            RawEdge *w_um = we + e_um;
            RawEdge *w_mv = we + e_mv;
            RawEdge *w_vm = we + e_vm;
            RawEdge *w_mu = we + e_mu;

            if (w_um->road_class != w_mv->road_class ||
                w_vm->road_class != w_mu->road_class ||
                (w_um->flags & EDGE_TOLL) != (w_mv->flags & EDGE_TOLL) ||
                (w_vm->flags & EDGE_TOLL) != (w_mu->flags & EDGE_TOLL)) {
                continue;
            }

            /* Merge u->m->v into e_um, kill e_mv */
            float len_uv = w_um->length_m + w_mv->length_m;
            w_um->f_speed  = (w_um->f_speed * w_um->length_m + w_mv->f_speed * w_mv->length_m) / len_uv;
            w_um->f_smooth = (w_um->f_smooth * w_um->length_m + w_mv->f_smooth * w_mv->length_m) / len_uv;
            w_um->length_m = len_uv;
            w_um->travel_time_s += w_mv->travel_time_s;
            w_um->to       = v;
            w_um->flags   |= EDGE_CONTRACTED | (w_mv->flags & EDGE_BRIDGE);
            *(in_next + e_um) = *(in_head + v);
            *(in_head + v)    = e_um;
            w_mv->alive = 0;

            /* Merge v->m->u into e_vm, kill e_mu */
            float len_vu = w_vm->length_m + w_mu->length_m;
            w_vm->f_speed  = (w_vm->f_speed * w_vm->length_m + w_mu->f_speed * w_mu->length_m) / len_vu;
            w_vm->f_smooth = (w_vm->f_smooth * w_vm->length_m + w_mu->f_smooth * w_mu->length_m) / len_vu;
            w_vm->length_m = len_vu;
            w_vm->travel_time_s += w_mu->travel_time_s;
            w_vm->to       = u;
            w_vm->flags   |= EDGE_CONTRACTED | (w_mu->flags & EDGE_BRIDGE);
            *(in_next + e_vm) = *(in_head + u);
            *(in_head + u)    = e_vm;
            w_mu->alive = 0;

            *(out_deg + m) = 0;
            *(in_deg + m)  = 0;
            *(removed + m) = 1;
            contracted_cnt++;
        }
    }
    xfree(cq);
    xfree(in_cq);
    xfree(out_head);
    xfree(in_head);
    xfree(out_next);
    xfree(in_next);
    xfree(out_deg);
    xfree(in_deg);
    xfree(no_contract);

    /* Step 3: Compact surviving nodes & live edges in place */
    uint32_t *remap = (uint32_t *)xmalloc((size_t)rn * sizeof(uint32_t));
    if (!remap) {
        xfree(removed);
        return -1;
    }
    uint32_t sn = 0;
    for (uint32_t i = 0; i < rn; i++) {
        if (!*(removed + i)) {
            *(remap + i) = sn++;
        } else {
            *(remap + i) = VMR_INVALID;
        }
    }
    xfree(removed);

    uint32_t sm = 0;
    for (uint32_t e = 0; e < rm; e++) {
        RawEdge *cur = we + e;
        if (!cur->alive) {
            continue;
        }
        uint32_t nu = *(remap + cur->from);
        uint32_t nv = *(remap + cur->to);
        if (nu == VMR_INVALID || nv == VMR_INVALID || nu == nv) {
            continue;
        }
        cur->from = nu;
        cur->to   = nv;
        if (sm != e) {
            we[sm] = *cur;
        }
        sm++;
    }

    if (sm > 1) {
        qsort(we, sm, sizeof(RawEdge), cmp_raw_edge);
    }

    g->node_count = sn;
    g->edge_count = sm;
    g->contracted_nodes = contracted_cnt;

    /* Allocate final SoA node & edge arrays */
    g->node_osm_id  = (int64_t *)xmalloc((size_t)sn * sizeof(int64_t));
    g->node_lat     = (double *)xmalloc((size_t)sn * sizeof(double));
    g->node_lon     = (double *)xmalloc((size_t)sn * sizeof(double));
    g->node_x       = (float *)xmalloc((size_t)sn * sizeof(float));
    g->node_y       = (float *)xmalloc((size_t)sn * sizeof(float));
    g->node_flags   = (uint8_t *)xcalloc(sn, sizeof(uint8_t));
    g->fwd_offsets  = (uint32_t *)xcalloc((size_t)sn + 1U, sizeof(uint32_t));

    size_t esz = (size_t)(sm > 0 ? sm : 1);
    g->edge_from       = (uint32_t *)xmalloc(esz * sizeof(uint32_t));
    g->edge_to         = (uint32_t *)xmalloc(esz * sizeof(uint32_t));
    g->edge_dx         = (float *)xmalloc(esz * sizeof(float));
    g->edge_dy         = (float *)xmalloc(esz * sizeof(float));
    g->edge_ux         = (float *)xmalloc(esz * sizeof(float));
    g->edge_uy         = (float *)xmalloc(esz * sizeof(float));
    g->edge_f_speed    = (float *)xmalloc(esz * sizeof(float));
    g->edge_f_class    = (float *)xmalloc(esz * sizeof(float));
    g->edge_f_straight = (float *)xmalloc(esz * sizeof(float));
    g->edge_f_smooth   = (float *)xmalloc(esz * sizeof(float));
    g->edge_length_m   = (float *)xmalloc(esz * sizeof(float));
    g->edge_time_s     = (float *)xmalloc(esz * sizeof(float));
    g->edge_flags      = (uint8_t *)xmalloc(esz * sizeof(uint8_t));

    if (!g->node_osm_id || !g->node_lat || !g->node_lon || !g->node_x || !g->node_y ||
        !g->node_flags || !g->fwd_offsets || !g->edge_from || !g->edge_to ||
        !g->edge_dx || !g->edge_dy || !g->edge_ux || !g->edge_uy ||
        !g->edge_f_speed || !g->edge_f_class || !g->edge_f_straight || !g->edge_f_smooth ||
        !g->edge_length_m || !g->edge_time_s || !g->edge_flags) {
        xfree(remap);
        return -1;
    }

    if (idmap_init(&g->osm_map, sn) != 0) {
        xfree(remap);
        return -1;
    }
    for (uint32_t i = 0; i < rn; i++) {
        uint32_t ni = *(remap + i);
        if (ni != VMR_INVALID) {
            int64_t oid = *(rg->node_osm_id + i);
            *(g->node_osm_id + ni) = oid;
            *(g->node_lat + ni)    = *(rg->node_lat + i);
            *(g->node_lon + ni)    = *(rg->node_lon + i);
            *(g->node_x + ni)      = *(rg->node_x + i);
            *(g->node_y + ni)      = *(rg->node_y + i);
            idmap_put(&g->osm_map, oid, ni);
        }
    }
    if (rg->turn_count > 0 && sn > 0) {
        uint32_t *cnt = (uint32_t *)xcalloc(sn, sizeof(uint32_t));
        uint32_t valid = 0;
        if (cnt) {
            for (uint32_t i = 0; i < rg->turn_count; i++) {
                uint32_t vi = *(remap + *(rg->turn_via + i));
                uint32_t fi = *(remap + *(rg->turn_from + i));
                uint32_t ti = *(remap + *(rg->turn_to + i));
                if (vi == VMR_INVALID || fi == VMR_INVALID || ti == VMR_INVALID || vi == fi || vi == ti) {
                    continue;
                }
                (*(cnt + vi))++;
                valid++;
            }
        }
        uint32_t *off = (uint32_t *)xcalloc((size_t)sn + 1U, sizeof(uint32_t));
        uint32_t *tf = NULL, *tt = NULL;
        uint8_t *tm = NULL;
        uint32_t *cur = NULL;
        if (off && valid > 0) {
            uint32_t sum = 0;
            for (uint32_t i = 0; i < sn; i++) {
                *(off + i) = sum;
                sum += *(cnt + i);
            }
            *(off + sn) = sum;
            tf = (uint32_t *)xmalloc((size_t)valid * sizeof(uint32_t));
            tt = (uint32_t *)xmalloc((size_t)valid * sizeof(uint32_t));
            tm = (uint8_t *)xmalloc((size_t)valid);
            cur = (uint32_t *)xmalloc((size_t)sn * sizeof(uint32_t));
            if (tf && tt && tm && cur) {
                memcpy(cur, off, (size_t)sn * sizeof(uint32_t));
                uint32_t w = 0;
                for (uint32_t i = 0; i < rg->turn_count; i++) {
                    uint32_t vi = *(remap + *(rg->turn_via + i));
                    uint32_t fi = *(remap + *(rg->turn_from + i));
                    uint32_t ti = *(remap + *(rg->turn_to + i));
                    if (vi == VMR_INVALID || fi == VMR_INVALID || ti == VMR_INVALID || vi == fi || vi == ti) {
                        continue;
                    }
                    uint32_t at = (*(cur + vi))++;
                    *(tf + at) = fi;
                    *(tt + at) = ti;
                    *(tm + at) = *(rg->turn_mode + i);
                    *(g->node_flags + vi) |= NODE_HAS_TURN_RES;
                    w++;
                }
                if (w == valid) {
                    g->turn_offsets = off;
                    g->turn_from = tf;
                    g->turn_to = tt;
                    g->turn_mode = tm;
                    g->turn_count = valid;
                    off = NULL;
                    tf = tt = NULL;
                    tm = NULL;
                }
            }
        }
        xfree(cnt);
        xfree(off);
        xfree(tf);
        xfree(tt);
        xfree(tm);
        xfree(cur);
    }
    xfree(remap);

    /* Raw node storage has been copied into the SoA graph: release it before
     * allocating/packing the much larger edge arrays. */
    xfree(rg->node_osm_id);
    xfree(rg->node_lat);
    xfree(rg->node_lon);
    xfree(rg->node_x);
    xfree(rg->node_y);
    idmap_free(&rg->osm_map);
    rg->node_osm_id = NULL;
    rg->node_lat = NULL;
    rg->node_lon = NULL;
    rg->node_x = NULL;
    rg->node_y = NULL;
    rg->node_count = 0;
    rg->node_cap = 0;

    /* Populate 8D SoA edge arrays & Forward CSR */
    for (uint32_t e = 0; e < sm; e++) {
        const RawEdge *w = we + e;
        uint32_t u = w->from;
        uint32_t v = w->to;

        *(g->edge_from + e) = u;
        *(g->edge_to + e)   = v;

        float dx = *(g->node_x + v) - *(g->node_x + u);
        float dy = *(g->node_y + v) - *(g->node_y + u);
        float mag = sqrtf(dx * dx + dy * dy);
        float ux = (mag > EPS_F) ? (dx / mag) : 0.0f;
        float uy = (mag > EPS_F) ? (dy / mag) : 0.0f;

        float straight = (w->length_m > EPS_F) ? (mag / w->length_m) : 1.0f;
        if (straight > 1.0f) straight = 1.0f;
        if (straight < 0.05f) straight = 0.05f;

        *(g->edge_dx + e)         = dx;
        *(g->edge_dy + e)         = dy;
        *(g->edge_ux + e)         = ux;
        *(g->edge_uy + e)         = uy;
        *(g->edge_f_speed + e)    = w->f_speed;
        *(g->edge_f_class + e)    = w->f_class;
        *(g->edge_f_straight + e) = straight;
        *(g->edge_f_smooth + e)   = w->f_smooth;
        *(g->edge_length_m + e)   = w->length_m;
        *(g->edge_time_s + e)     = w->travel_time_s;
        *(g->edge_flags + e)      = w->flags;

        (*(g->fwd_offsets + u + 1U))++;
    }

    /* Precompute attr_scale for the five built-in profiles: removes the 4D dot
     * product and division from the query hot loop. */
    static const float pw[5][4] = {
        {1.20f, 1.00f, 0.60f, 0.50f},
        {2.20f, 1.50f, 0.30f, 0.40f},
        {0.50f, 0.30f, 1.20f, 0.30f},
        {1.80f, 2.20f, 0.50f, 0.60f},
        {0.70f, 0.20f, 0.80f, 1.40f}
    };
    if (sm > 0) {
        for (uint32_t p = 0; p < 5; p++) {
            float *arr = (float *)xmalloc((size_t)sm * sizeof(float));
            if (!arr) {
                continue;
            }
            for (uint32_t e = 0; e < sm; e++) {
                float dot = pw[p][0] * *(g->edge_f_speed + e) +
                            pw[p][1] * *(g->edge_f_class + e) +
                            pw[p][2] * *(g->edge_f_straight + e) +
                            pw[p][3] * *(g->edge_f_smooth + e);
                *(arr + e) = 1.0f / (1.0f + 0.22f * dot);
            }
            g->edge_attr_scale[p] = arr;
        }
    }

    /* Raw edge pool is now mirrored in the SoA arrays. */
    xfree(rg->edges);
    rg->edges = NULL;
    rg->edge_count = 0;
    rg->edge_cap = 0;

    for (uint32_t i = 1; i <= sn; i++) {
        *(g->fwd_offsets + i) += *(g->fwd_offsets + i - 1U);
    }

    if (vmr_build_reverse_csr_and_components(g) != 0) {
        return -1;
    }

    /* Step 4: Dead-end 1-core cul-de-sac peeling.
     * Degree is always the number of DISTINCT undirected neighbours: a two-way
     * edge counts once, no matter whether the node is reached through fwd or
     * rev CSR. The stamp array guarantees each neighbour is decremented at most
     * once per peeled node (previously a two-way stub decremented a junction
     * twice and could falsely mark a live junction as a dead end). */
    uint32_t *undir_deg  = (uint32_t *)xcalloc(sn, sizeof(uint32_t));
    uint32_t *nbr_stamp  = (uint32_t *)xcalloc(sn, sizeof(uint32_t));
    uint32_t *peel_q     = (uint32_t *)xmalloc((size_t)sn * sizeof(uint32_t));
    if (undir_deg && nbr_stamp && peel_q) {
        uint32_t pq_h = 0, pq_t = 0;
        for (uint32_t u = 0; u < sn; u++) {
            uint32_t distinct = 0;
            for (uint32_t ei = *(g->fwd_offsets + u); ei < *(g->fwd_offsets + u + 1U); ei++) {
                uint32_t v = *(g->edge_to + ei);
                if (*(nbr_stamp + v) != u + 1U) {
                    *(nbr_stamp + v) = u + 1U;
                    distinct++;
                }
            }
            for (uint32_t ri = *(g->rev_offsets + u); ri < *(g->rev_offsets + u + 1U); ri++) {
                uint32_t v = *(g->rev_src_node + ri);
                if (*(nbr_stamp + v) != u + 1U) {
                    *(nbr_stamp + v) = u + 1U;
                    distinct++;
                }
            }
            *(undir_deg + u) = distinct;
            if (distinct <= 1) {
                *(g->node_flags + u) |= NODE_DEAD_TREE;
                *(peel_q + pq_t++) = u;
            } else {
                *(g->node_flags + u) |= NODE_JUNCTION;
            }
        }
        /* Reuse the stamp array for per-peel neighbour de-duplication. */
        memset(nbr_stamp, 0, (size_t)sn * sizeof(uint32_t));

        uint32_t peeled = 0;
        while (pq_h < pq_t) {
            uint32_t u = *(peel_q + pq_h++);
            uint32_t ustamp = u + 1U;
            peeled++;
            for (uint32_t ei = *(g->fwd_offsets + u); ei < *(g->fwd_offsets + u + 1U); ei++) {
                uint32_t v = *(g->edge_to + ei);
                if (*(nbr_stamp + v) == ustamp) continue;
                *(nbr_stamp + v) = ustamp;
                if ((*(g->node_flags + v) & NODE_DEAD_TREE) || *(undir_deg + v) == 0) continue;
                if (--(*(undir_deg + v)) <= 1) {
                    *(g->node_flags + v) |= NODE_DEAD_TREE;
                    *(peel_q + pq_t++) = v;
                }
            }
            for (uint32_t ri = *(g->rev_offsets + u); ri < *(g->rev_offsets + u + 1U); ri++) {
                uint32_t v = *(g->rev_src_node + ri);
                if (*(nbr_stamp + v) == ustamp) continue;
                *(nbr_stamp + v) = ustamp;
                if ((*(g->node_flags + v) & NODE_DEAD_TREE) || *(undir_deg + v) == 0) continue;
                if (--(*(undir_deg + v)) <= 1) {
                    *(g->node_flags + v) |= NODE_DEAD_TREE;
                    *(peel_q + pq_t++) = v;
                }
            }
        }
        g->peeled_dead_nodes = peeled;
    }
    xfree(undir_deg);
    xfree(nbr_stamp);
    xfree(peel_q);

    for (uint32_t e = 0; e < g->edge_count; e++) {
        if (*(g->node_flags + *(g->edge_from + e)) & NODE_DEAD_TREE) {
            *(g->edge_flags + e) |= EDGE_DEAD_SOURCE;
        }
        if (*(g->node_flags + *(g->edge_to + e)) & NODE_DEAD_TREE) {
            *(g->edge_flags + e) |= EDGE_DEAD_TARGET;
        }
    }

    vmr_reorder_graph(g);

    vmr_build_two_hop(g);
    vmr_build_two_hop_rev(g);
    if (vmr_rev_enabled()) {
        vmr_build_edge_blocks(g);
    }

    if (vmr_build_spatial_grid(g) != 0) {
        return -1;
    }

    g->landmark_count = 32;
    vmr_build_landmarks(g);

    /* Contraction Hierarchy (CCH) on travel time; opt-in via VMR_CCH so the
     * default vector pipeline stays lean. */
    if (getenv("VMR_CCH") != NULL) {
        clock_t tc0 = clock();
        if (vmr_build_ch(g) != 0) {
            LOG_WARN("CCH build failed; CCH mode disabled");
        } else {
            LOG_INFO("CCH built in %.1f ms (%u levels, %u ch-edges)",
                     1000.0 * (double)(clock() - tc0) / (double)CLOCKS_PER_SEC,
                     g->node_count, g->ch_edge_count);
        }
    }

    /* Vector CH is not used for the exact-vector A* (UK measured slower). */

    clock_t t1 = clock();
    double ms = 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
    vmr_compute_attr_min(g);
    LOG_INFO("Preprocessed UK SoA Graph in %.1f ms: %u -> %u nodes (-%.1f%% contracted), "
             "%u -> %u edges, %u cul-de-sac nodes peeled",
             ms, rn, sn,
             100.0 * (double)contracted_cnt / (double)rn,
             rm, sm, g->peeled_dead_nodes);
    return 0;
}

/* ========================================================================= */
/* SECTION 6: BINARY .vgraph CACHE SERIALIZATION                              */
/* ========================================================================= */

typedef struct {
    uint32_t magic;
    uint32_t version;
    uint32_t node_count;
    uint32_t edge_count;
    uint32_t raw_node_count;
    uint32_t raw_edge_count;
    uint32_t contracted_nodes;
    uint32_t peeled_dead_nodes;
    uint32_t grid_dim;
    float    grid_min_x;
    float    grid_max_x;
    float    grid_min_y;
    float    grid_max_y;
} VgraphHeader;

static int vmr_fwrite_chunk(const void *ptr, size_t sz, size_t n, FILE *fp) {
    size_t off = 0;
    while (off < n) {
        size_t cnt = n - off;
        if (cnt > 262144U) {
            cnt = 262144U;
        }
        if (fwrite((const uint8_t *)ptr + off * sz, sz, cnt, fp) != cnt) {
            return 0;
        }
        off += cnt;
    }
    return 1;
}

static int vmr_save_binary(const VmrGraph *g, const char *path) {
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        LOG_ERR("Cannot open '%s' for writing", path);
        return -1;
    }
    VgraphHeader hdr;
    memset(&hdr, 0, sizeof(hdr));
    hdr.magic = VGRAPH_MAGIC;
    hdr.version = VGRAPH_FORMAT_V;
    hdr.node_count = g->node_count;
    hdr.edge_count = g->edge_count;
    hdr.raw_node_count = g->raw_node_count;
    hdr.raw_edge_count = g->raw_edge_count;
    hdr.contracted_nodes = g->contracted_nodes;
    hdr.peeled_dead_nodes = g->peeled_dead_nodes;
    hdr.grid_dim = SPATIAL_GRID_DIM;
    hdr.grid_min_x = g->grid_min_x;
    hdr.grid_max_x = g->grid_max_x;
    hdr.grid_min_y = g->grid_min_y;
    hdr.grid_max_y = g->grid_max_y;

    size_t n = g->node_count;
    size_t m = g->edge_count;
    size_t cells = (size_t)(SPATIAL_GRID_DIM * SPATIAL_GRID_DIM);

    int ok = 1;
    ok &= (fwrite(&hdr, sizeof(hdr), 1, fp) == 1);
    ok &= (fwrite(g->node_osm_id,  sizeof(int64_t),  n, fp) == n);
    ok &= (fwrite(g->node_lat,     sizeof(double),   n, fp) == n);
    ok &= (fwrite(g->node_lon,     sizeof(double),   n, fp) == n);
    ok &= (fwrite(g->node_x,       sizeof(float),    n, fp) == n);
    ok &= (fwrite(g->node_y,       sizeof(float),    n, fp) == n);
    ok &= (fwrite(g->node_comp_id, sizeof(uint32_t), n, fp) == n);
    ok &= (fwrite(g->node_flags,   sizeof(uint8_t),  n, fp) == n);
    ok &= (fwrite(g->fwd_offsets,  sizeof(uint32_t), n + 1U, fp) == n + 1U);

    if (m > 0) {
        ok &= (fwrite(g->edge_from,       sizeof(uint32_t), m, fp) == m);
        ok &= (fwrite(g->edge_to,         sizeof(uint32_t), m, fp) == m);
        ok &= (fwrite(g->edge_dx,         sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_dy,         sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_ux,         sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_uy,         sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_f_speed,    sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_f_class,    sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_f_straight, sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_f_smooth,   sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_length_m,   sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_time_s,     sizeof(float),    m, fp) == m);
        ok &= (fwrite(g->edge_flags,      sizeof(uint8_t),  m, fp) == m);
    }

    uint8_t amask = 0;
    for (uint32_t p = 0; p < 5; p++) {
        if (g->edge_attr_scale[p]) {
            amask |= (uint8_t)(1U << p);
        }
    }
    ok &= (fwrite(&amask, sizeof(uint8_t), 1, fp) == 1);
    for (uint32_t p = 0; p < 5; p++) {
        if (amask & (1U << p)) {
            ok &= (fwrite(g->edge_attr_scale[p], sizeof(float), m, fp) == m);
        }
    }

    ok &= (fwrite(g->rev_offsets,  sizeof(uint32_t), n + 1U, fp) == n + 1U);
    if (m > 0) {
        ok &= (fwrite(g->rev_src_node, sizeof(uint32_t), m, fp) == m);
        ok &= (fwrite(g->rev_fwd_eid,  sizeof(uint32_t), m, fp) == m);
    }
    ok &= (fwrite(g->grid_offsets, sizeof(uint32_t), cells + 1U, fp) == cells + 1U);
    ok &= (fwrite(g->grid_nodes,   sizeof(uint32_t), n, fp) == n);

    uint32_t l = g->landmark_count;
    ok &= (fwrite(&l, sizeof(uint32_t), 1, fp) == 1);
    if (l > 0 && n > 0) {
        ok &= (fwrite(g->landmark_node, sizeof(uint32_t), l, fp) == l);
        ok &= vmr_fwrite_chunk(g->alt_nodes, sizeof(AltNode32), n, fp);
    }

    uint32_t tc = g->turn_count;
    ok &= (fwrite(&tc, sizeof(uint32_t), 1, fp) == 1);
    if (tc > 0 && n > 0) {
        ok &= (fwrite(g->turn_offsets, sizeof(uint32_t), n + 1U, fp) == n + 1U);
        ok &= (fwrite(g->turn_from,    sizeof(uint32_t), tc, fp) == tc);
        ok &= (fwrite(g->turn_to,      sizeof(uint32_t), tc, fp) == tc);
        ok &= (fwrite(g->turn_mode,    sizeof(uint8_t),  tc, fp) == tc);
    }

    uint32_t cb = g->cch_built ? 1U : 0U;
    ok &= (fwrite(&cb, sizeof(uint32_t), 1, fp) == 1);
    if (cb && n > 0) {
        uint32_t ce = g->ch_edge_count;
        ok &= (fwrite(&ce, sizeof(uint32_t), 1, fp) == 1);
        ok &= (fwrite(g->ch_level, sizeof(uint32_t), n, fp) == n);
        ok &= (fwrite(g->ch_from, sizeof(uint32_t), ce, fp) == ce);
        ok &= (fwrite(g->ch_to, sizeof(uint32_t), ce, fp) == ce);
        ok &= (fwrite(g->ch_time, sizeof(float), ce, fp) == ce);
        ok &= (fwrite(g->ch_fwd_off, sizeof(uint32_t), n + 1U, fp) == n + 1U);
        ok &= (fwrite(g->ch_rev_off, sizeof(uint32_t), n + 1U, fp) == n + 1U);
        ok &= (fwrite(g->ch_rev_src, sizeof(uint32_t), ce, fp) == ce);
        ok &= (fwrite(g->ch_rev_eid, sizeof(uint32_t), ce, fp) == ce);
        ok &= (fwrite(g->ch_via, sizeof(uint32_t), ce, fp) == ce);
        ok &= (fwrite(g->ch_child1, sizeof(uint32_t), ce, fp) == ce);
        ok &= (fwrite(g->ch_child2, sizeof(uint32_t), ce, fp) == ce);
    }

    /* Persist the two-hop continuation tables so a reload can mmap them
     * instead of rebuilding (saves ~360 MB of heap and the O(E) rebuild). */
    uint32_t thb = (g->th_off != NULL && m > 0) ? 1U : 0U;
    ok &= (fwrite(&thb, sizeof(uint32_t), 1, fp) == 1);
    if (thb) {
        uint32_t thc = g->th_off[m];
        ok &= (fwrite(&thc, sizeof(uint32_t), 1, fp) == 1);
        ok &= (fwrite(g->th_off, sizeof(uint32_t), m + 1U, fp) == m + 1U);
        ok &= (fwrite(g->th_sx, sizeof(float), thc, fp) == thc);
        ok &= (fwrite(g->th_sy, sizeof(float), thc, fp) == thc);
        ok &= (fwrite(g->th_inv, sizeof(float), thc, fp) == thc);
    }
    uint32_t rthb = (g->rth_off != NULL && m > 0) ? 1U : 0U;
    ok &= (fwrite(&rthb, sizeof(uint32_t), 1, fp) == 1);
    if (rthb) {
        uint32_t rthc = g->rth_off[m];
        ok &= (fwrite(&rthc, sizeof(uint32_t), 1, fp) == 1);
        ok &= (fwrite(g->rth_off, sizeof(uint32_t), m + 1U, fp) == m + 1U);
        ok &= (fwrite(g->rth_sx, sizeof(float), rthc, fp) == rthc);
        ok &= (fwrite(g->rth_sy, sizeof(float), rthc, fp) == rthc);
        ok &= (fwrite(g->rth_inv, sizeof(float), rthc, fp) == rthc);
    }

    fclose(fp);
    if (!ok) {
        LOG_ERR("I/O error writing '%s'", path);
        return -1;
    }
    LOG_INFO("Saved preprocessed binary SoA graph to '%s' (%u nodes, %u edges)",
             path, g->node_count, g->edge_count);
    return 0;
}

static int vmr_load_binary(VmrGraph *g, const char *path) {
    vmr_graph_free(g);
    vmr_graph_init(g);

    FILE *fp = fopen(path, "rb");
    if (fp == NULL) {
        LOG_ERR("Cannot open '%s' for reading", path);
        return -1;
    }
    VgraphHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, fp) != 1 ||
        hdr.magic != VGRAPH_MAGIC || hdr.version != VGRAPH_FORMAT_V ||
        hdr.grid_dim != SPATIAL_GRID_DIM) {
        LOG_ERR("Invalid or incompatible .vgraph header in '%s'", path);
        fclose(fp);
        return -1;
    }

    clock_t t0 = clock();
    size_t n = hdr.node_count;
    size_t m = hdr.edge_count;
    size_t esz = (m > 0) ? m : 1;
    size_t cells = (size_t)(SPATIAL_GRID_DIM * SPATIAL_GRID_DIM);

    g->node_count = hdr.node_count;
    g->edge_count = hdr.edge_count;
    g->raw_node_count = hdr.raw_node_count;
    g->raw_edge_count = hdr.raw_edge_count;
    g->contracted_nodes = hdr.contracted_nodes;
    g->peeled_dead_nodes = hdr.peeled_dead_nodes;
    g->grid_min_x = hdr.grid_min_x;
    g->grid_max_x = hdr.grid_max_x;
    g->grid_min_y = hdr.grid_min_y;
    g->grid_max_y = hdr.grid_max_y;

    {
        HANDLE fh = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ, NULL,
                                OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
        if (fh == INVALID_HANDLE_VALUE) {
            fclose(fp);
            vmr_graph_free(g);
            return -1;
        }
        LARGE_INTEGER fsz;
        if (!GetFileSizeEx(fh, &fsz) || fsz.QuadPart < (LONGLONG)sizeof(VgraphHeader)) {
            CloseHandle(fh);
            fclose(fp);
            vmr_graph_free(g);
            return -1;
        }
        HANDLE mh = CreateFileMappingA(fh, NULL, PAGE_READONLY, 0, 0, NULL);
        CloseHandle(fh);
        if (!mh) {
            fclose(fp);
            vmr_graph_free(g);
            return -1;
        }
        void *base = MapViewOfFile(mh, FILE_MAP_READ, 0, 0, 0);
        if (!base) {
            CloseHandle(mh);
            fclose(fp);
            vmr_graph_free(g);
            return -1;
        }
        g->map = base;
        g->map_h = mh;
        uint8_t *p = (uint8_t *)base;
        if (memcmp(p, &hdr, sizeof(hdr)) != 0) {
            /* header was already read; mapped copy must match */
        }
        size_t off = sizeof(VgraphHeader);
        g->node_osm_id  = (int64_t *)(p + off); off += (size_t)n * sizeof(int64_t);
        g->node_lat     = (double  *)(p + off); off += (size_t)n * sizeof(double);
        g->node_lon     = (double  *)(p + off); off += (size_t)n * sizeof(double);
        g->node_x       = (float   *)(p + off); off += (size_t)n * sizeof(float);
        g->node_y       = (float   *)(p + off); off += (size_t)n * sizeof(float);
        g->node_comp_id = (uint32_t*)(p + off); off += (size_t)n * sizeof(uint32_t);
        g->node_flags   = (uint8_t *)(p + off); off += (size_t)n * sizeof(uint8_t);
        g->fwd_offsets  = (uint32_t*)(p + off); off += (size_t)(n + 1U) * sizeof(uint32_t);
        g->edge_from       = (uint32_t*)(p + off); off += esz * sizeof(uint32_t);
        g->edge_to         = (uint32_t*)(p + off); off += esz * sizeof(uint32_t);
        g->edge_dx         = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_dy         = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_ux         = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_uy         = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_f_speed    = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_f_class    = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_f_straight = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_f_smooth   = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_length_m   = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_time_s     = (float   *)(p + off); off += esz * sizeof(float);
        g->edge_flags      = (uint8_t *)(p + off); off += esz * sizeof(uint8_t);

        /* attr_scale block (mask + float[m] per set profile) sits between the
         * edge arrays and the reverse CSR on disk; consume it so the mapped
         * suffix stays aligned. Copied to owned buffers (freed unconditionally). */
        uint8_t amask = 0;
        if (off < (size_t)fsz.QuadPart) { amask = *(const uint8_t *)(p + off); off += 1U; }
        for (uint32_t pi = 0; pi < 5; pi++) {
            if (amask & (uint8_t)(1U << pi)) {
                g->edge_attr_scale[pi] = (float *)xmalloc((size_t)(m ? m : 1) * sizeof(float));
                if (g->edge_attr_scale[pi] && m > 0 &&
                    off + (size_t)m * sizeof(float) <= (size_t)fsz.QuadPart) {
                    memcpy(g->edge_attr_scale[pi], p + off, (size_t)m * sizeof(float));
                }
                off += (size_t)m * sizeof(float);
            }
        }

        g->rev_offsets  = (uint32_t*)(p + off); off += (size_t)(n + 1U) * sizeof(uint32_t);
        g->rev_src_node = (uint32_t*)(p + off); off += esz * sizeof(uint32_t);
        g->rev_fwd_eid  = (uint32_t*)(p + off); off += esz * sizeof(uint32_t);
        g->grid_offsets = (uint32_t*)(p + off); off += (size_t)(cells + 1U) * sizeof(uint32_t);
        g->grid_nodes   = (uint32_t*)(p + off); off += (size_t)n * sizeof(uint32_t);

        uint32_t l = 0;
        if (off + 4 <= (size_t)fsz.QuadPart) { l = *(const uint32_t *)(p + off); off += 4U; }
        if (l > 64U) { l = 0; }
        g->landmark_count = l;
        if (l > 0 && n > 0) {
            /* point both into the mapping (no heap, no 124-byte leak) */
            g->landmark_node = (uint32_t *)(p + off);
            off += (size_t)l * sizeof(uint32_t);
            /* mmap the (~695 MB) ALT table in place instead of a heap copy */
            g->alt_nodes = (AltNode32 *)(p + off);
            off += (size_t)n * sizeof(AltNode32);
        }

        uint32_t tc = 0;
        if (off + 4 <= (size_t)fsz.QuadPart) { tc = *(const uint32_t *)(p + off); off += 4U; }
        if (tc > 100000000U) { tc = 0; }
        g->turn_count = tc;
        if (tc > 0 && n > 0) {
            g->turn_offsets = (uint32_t *)xcalloc((size_t)n + 1U, sizeof(uint32_t));
            g->turn_from = (uint32_t *)xmalloc((size_t)tc * sizeof(uint32_t));
            g->turn_to   = (uint32_t *)xmalloc((size_t)tc * sizeof(uint32_t));
            g->turn_mode = (uint8_t  *)xmalloc((size_t)tc);
            if (g->turn_offsets && g->turn_from && g->turn_to && g->turn_mode &&
                off + (size_t)(n + 1U) * sizeof(uint32_t) + (size_t)tc * sizeof(uint32_t) * 2U + (size_t)tc <= (size_t)fsz.QuadPart) {
                memcpy(g->turn_offsets, p + off, (size_t)(n + 1U) * sizeof(uint32_t)); off += (size_t)(n + 1U) * sizeof(uint32_t);
                memcpy(g->turn_from, p + off, (size_t)tc * sizeof(uint32_t)); off += (size_t)tc * sizeof(uint32_t);
                memcpy(g->turn_to,   p + off, (size_t)tc * sizeof(uint32_t)); off += (size_t)tc * sizeof(uint32_t);
                memcpy(g->turn_mode, p + off, (size_t)tc); off += (size_t)tc;
            }
        }

        uint32_t cb = 0;
        if (off + 4 <= (size_t)fsz.QuadPart) { cb = *(const uint32_t *)(p + off); off += 4U; }
        if (cb && n > 0 && off + 4 <= (size_t)fsz.QuadPart) {
            uint32_t ce = *(const uint32_t *)(p + off); off += 4U;
            size_t need = (size_t)n * 4 + ((size_t)ce * 4) * 3 + ((size_t)n + 1U) * 4 * 2 + ((size_t)ce * 4) * 5;
            if (off + need <= (size_t)fsz.QuadPart) {
                g->ch_edge_count = ce;
                /* mmap the CH arrays in place (owned by the mapping, paged on demand) */
                g->ch_level   = (uint32_t *)(p + off); off += (size_t)n * sizeof(uint32_t);
                g->ch_from    = (uint32_t *)(p + off); off += (size_t)ce * sizeof(uint32_t);
                g->ch_to      = (uint32_t *)(p + off); off += (size_t)ce * sizeof(uint32_t);
                g->ch_time    = (float   *)(p + off); off += (size_t)ce * sizeof(float);
                g->ch_fwd_off = (uint32_t *)(p + off); off += ((size_t)n + 1U) * sizeof(uint32_t);
                g->ch_rev_off = (uint32_t *)(p + off); off += ((size_t)n + 1U) * sizeof(uint32_t);
                g->ch_rev_src = (uint32_t *)(p + off); off += (size_t)ce * sizeof(uint32_t);
                g->ch_rev_eid = (uint32_t *)(p + off); off += (size_t)ce * sizeof(uint32_t);
                g->ch_via     = (uint32_t *)(p + off); off += (size_t)ce * sizeof(uint32_t);
                g->ch_child1  = (uint32_t *)(p + off); off += (size_t)ce * sizeof(uint32_t);
                g->ch_child2  = (uint32_t *)(p + off); off += (size_t)ce * sizeof(uint32_t);
                g->cch_built = 1U;
            }
        }

        /* Persisted two-hop continuation tables: mmap in place (no heap, no
         * rebuild). Absent in older caches -> fall back to building them. */
        uint32_t thb = 0;
        if (off + 4 <= (size_t)fsz.QuadPart) { thb = *(const uint32_t *)(p + off); off += 4U; }
        if (thb && off + 4 <= (size_t)fsz.QuadPart) {
            uint32_t thc = *(const uint32_t *)(p + off); off += 4U;
            g->th_off = (uint32_t *)(p + off); off += (size_t)(m + 1U) * sizeof(uint32_t);
            g->th_sx = (float *)(p + off); off += (size_t)thc * sizeof(float);
            g->th_sy = (float *)(p + off); off += (size_t)thc * sizeof(float);
            g->th_inv = (float *)(p + off); off += (size_t)thc * sizeof(float);
            g->th_count = thc;
        }
        uint32_t rthb = 0;
        if (off + 4 <= (size_t)fsz.QuadPart) { rthb = *(const uint32_t *)(p + off); off += 4U; }
        if (rthb && off + 4 <= (size_t)fsz.QuadPart) {
            uint32_t rthc = *(const uint32_t *)(p + off); off += 4U;
            g->rth_off = (uint32_t *)(p + off); off += (size_t)(m + 1U) * sizeof(uint32_t);
            g->rth_sx = (float *)(p + off); off += (size_t)rthc * sizeof(float);
            g->rth_sy = (float *)(p + off); off += (size_t)rthc * sizeof(float);
            g->rth_inv = (float *)(p + off); off += (size_t)rthc * sizeof(float);
            g->rth_count = rthc;
        }
        fclose(fp);
    }

    int ok = 1;
    (void)ok;

    if (!ok || idmap_init(&g->osm_map, n) != 0) {
        vmr_graph_free(g);
        return -1;
    }
    const int64_t *oids = g->node_osm_id;
    for (uint32_t i = 0; i < (uint32_t)n; i++) {
        idmap_put(&g->osm_map, *(oids + i), i);
    }

    if (!g->th_off)  vmr_build_two_hop(g);
    if (!g->rth_off) vmr_build_two_hop_rev(g);
    if (vmr_rev_enabled()) {
        vmr_build_edge_blocks(g);
    }

    clock_t t1 = clock();
    vmr_compute_attr_min(g);
    LOG_INFO("Loaded binary SoA graph '%s' in %.1f ms (%zu nodes, %zu edges)",
             path, 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC, n, m);
    return 0;
}

/* ========================================================================= */
/* SECTION 7: SQLITE3 STORAGE, OSM XML UK EXTRACTOR & UK DEMO SEEDER         */
/* ========================================================================= */

static int sql_exec(sqlite3 *db, const char *sql) {
    char *err = NULL;
    if (sqlite3_exec(db, sql, NULL, NULL, &err) != SQLITE_OK) {
        LOG_ERR("SQLite error: %s", err ? err : sqlite3_errmsg(db));
        if (err) sqlite3_free(err);
        return -1;
    }
    return 0;
}

static int db_init_schema(sqlite3 *db) {
    const char *sql =
        "PRAGMA journal_mode = WAL;\n"
        "CREATE TABLE IF NOT EXISTS osm_nodes (\n"
        "    id  INTEGER PRIMARY KEY,\n"
        "    lat REAL NOT NULL,\n"
        "    lon REAL NOT NULL\n"
        ");\n"
        "CREATE TABLE IF NOT EXISTS osm_edges (\n"
        "    from_node  INTEGER NOT NULL,\n"
        "    to_node    INTEGER NOT NULL,\n"
        "    road_class INTEGER NOT NULL,\n"
        "    speed_kmh  REAL NOT NULL,\n"
        "    smooth     REAL NOT NULL DEFAULT 0.85,\n"
        "    flags      INTEGER NOT NULL DEFAULT 0\n"
        ");\n"
        "CREATE TABLE IF NOT EXISTS turn_restrictions (\n"
        "    via_node  INTEGER NOT NULL,\n"
        "    from_node INTEGER NOT NULL,\n"
        "    to_node   INTEGER NOT NULL,\n"
        "    mode      INTEGER NOT NULL DEFAULT 0,\n"
        "    PRIMARY KEY(via_node, from_node, to_node, mode)\n"
        ");\n";
    return sql_exec(db, sql);
}

static void road_class_defaults(RoadClass rc, float *speed_kmh, float *f_class) {
    switch (rc) {
        case RC_MOTORWAY:     *speed_kmh = 112.6f; *f_class = 1.00f; break; /* 70 mph */
        case RC_TRUNK:        *speed_kmh =  96.5f; *f_class = 0.85f; break; /* 60 mph */
        case RC_PRIMARY:      *speed_kmh =  80.5f; *f_class = 0.70f; break; /* 50 mph */
        case RC_SECONDARY:    *speed_kmh =  64.4f; *f_class = 0.55f; break; /* 40 mph */
        case RC_TERTIARY:     *speed_kmh =  48.3f; *f_class = 0.40f; break; /* 30 mph */
        case RC_RESIDENTIAL:  *speed_kmh =  32.2f; *f_class = 0.25f; break; /* 20 mph */
        default:              *speed_kmh =  48.3f; *f_class = 0.30f; break;
    }
}

static int parse_uk_highway_tag(const char *val, RoadClass *out_rc) {
    if (strcmp(val, "motorway") == 0 || strcmp(val, "motorway_link") == 0) {
        *out_rc = RC_MOTORWAY; return 1;
    }
    if (strcmp(val, "trunk") == 0 || strcmp(val, "trunk_link") == 0) {
        *out_rc = RC_TRUNK; return 1;
    }
    if (strcmp(val, "primary") == 0 || strcmp(val, "primary_link") == 0) {
        *out_rc = RC_PRIMARY; return 1;
    }
    if (strcmp(val, "secondary") == 0 || strcmp(val, "secondary_link") == 0) {
        *out_rc = RC_SECONDARY; return 1;
    }
    if (strcmp(val, "tertiary") == 0 || strcmp(val, "tertiary_link") == 0) {
        *out_rc = RC_TERTIARY; return 1;
    }
    if (strcmp(val, "residential") == 0 || strcmp(val, "living_street") == 0) {
        *out_rc = RC_RESIDENTIAL; return 1;
    }
    if (strcmp(val, "unclassified") == 0) {
        *out_rc = RC_UNCLASSIFIED; return 1;
    }
    return 0;
}

static float parse_uk_maxspeed_kmh(const char *val, RoadClass rc) {
    float def_spd = 48.3f, dummy = 0.0f;
    road_class_defaults(rc, &def_spd, &dummy);

    if (val == NULL || *val == '\0') {
        return def_spd;
    }
    if (c_stricmp(val, "national") == 0 || c_stricmp(val, "NSL") == 0) {
        return (rc == RC_MOTORWAY) ? 112.65f : 96.56f;
    }
    char *endp = NULL;
    double num = strtod(val, &endp);
    if (endp == val || num <= 0.0) {
        return def_spd;
    }
    while (*endp != '\0' && isspace((unsigned char)*endp)) {
        endp++;
    }
    if (c_stricmp(endp, "mph") == 0) {
        num *= 1.609344;
    }
    if (num < 8.0) num = 8.0;
    if (num > 140.0) num = 140.0;
    return (float)num;
}

/* Tiny XML attribute extractor: finds `attr="..."` or `attr='...'` on a line */
static int xml_get_attr(const char *line, const char *attr, char *out, size_t out_sz) {
    size_t alen = strlen(attr);
    const char *p = line;
    while ((p = strstr(p, attr)) != NULL) {
        if ((p == line || isspace((unsigned char)*(p - 1))) && *(p + alen) == '=') {
            const char *q = p + alen + 1;
            char quote = *q;
            if (quote == '"' || quote == '\'') {
                q++;
                const char *end = strchr(q, quote);
                if (end != NULL) {
                    size_t len = (size_t)(end - q);
                    if (len >= out_sz) len = out_sz - 1;
                    memcpy(out, q, len);
                    *(out + len) = '\0';
                    return 1;
                }
            }
        }
        p += alen;
    }
    return 0;
}

/*
 * Two-pass streaming OSM XML (.osm) importer tailored for Geofabrik UK extracts:
 *   Pass 1: Reads <way> elements, identifies routable highway=* ways, and marks
 *           referenced node IDs in a compact hash set so non-road nodes
 *           (buildings, trees, coastline) are discarded.
 *   Pass 2: Reads <node> elements for marked IDs and <way> elements to bulk
 *           insert edges into SQLite inside fast transactions.
 */
static int db_import_osm_xml(sqlite3 *db, const char *osm_path) {
    FILE *fp = fopen(osm_path, "r");
    if (fp == NULL) {
        LOG_ERR("Cannot open OSM XML file '%s'", osm_path);
        return -1;
    }

    IdMap needed_nodes;
    if (idmap_init(&needed_nodes, 65536) != 0) {
        fclose(fp);
        return -1;
    }

    char line[4096];
    int64_t *way_refs = (int64_t *)xmalloc(4096 * sizeof(int64_t));
    size_t way_ref_cap = 4096;
    size_t way_ref_cnt = 0;
    int in_way = 0;
    int is_highway = 0;
    size_t pass1_ways = 0;

    LOG_INFO("OSM Import Pass 1: scanning highway ways in '%s'...", osm_path);
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (!in_way) {
            if (strstr(line, "<way ") != NULL) {
                in_way = 1;
                is_highway = 0;
                way_ref_cnt = 0;
                if (strstr(line, "/>") != NULL) {
                    in_way = 0;
                }
            }
        } else {
            char val[128], key[128];
            if (strstr(line, "<nd ") != NULL && xml_get_attr(line, "ref", val, sizeof(val))) {
                if (way_ref_cnt >= way_ref_cap) {
                    way_ref_cap *= 2;
                    int64_t *nr = (int64_t *)xrealloc(way_refs, way_ref_cap * sizeof(int64_t));
                    if (!nr) break;
                    way_refs = nr;
                }
                *(way_refs + way_ref_cnt++) = strtoll(val, NULL, 10);
            } else if (strstr(line, "<tag ") != NULL &&
                       xml_get_attr(line, "k", key, sizeof(key)) &&
                       xml_get_attr(line, "v", val, sizeof(val))) {
                if (strcmp(key, "highway") == 0) {
                    RoadClass dummy_rc;
                    if (parse_uk_highway_tag(val, &dummy_rc)) {
                        is_highway = 1;
                    }
                }
            } else if (strstr(line, "</way>") != NULL) {
                if (is_highway && way_ref_cnt >= 2) {
                    pass1_ways++;
                    for (size_t i = 0; i < way_ref_cnt; i++) {
                        idmap_put(&needed_nodes, *(way_refs + i), 1U);
                    }
                }
                in_way = 0;
            }
        }
    }
    LOG_INFO("OSM Pass 1 complete: %zu routable ways, %zu unique highway nodes",
             pass1_ways, needed_nodes.count);

    rewind(fp);
    sql_exec(db, "PRAGMA synchronous = OFF;");
    sql_exec(db, "BEGIN TRANSACTION;");

    sqlite3_stmt *nstmt = NULL, *estmt = NULL;
    if (sqlite3_prepare_v2(db, "INSERT OR IGNORE INTO osm_nodes(id, lat, lon) VALUES(?,?,?);",
                           -1, &nstmt, NULL) != SQLITE_OK ||
        sqlite3_prepare_v2(db,
            "INSERT INTO osm_edges(from_node, to_node, road_class, speed_kmh, smooth, flags) "
            "VALUES(?,?,?,?,?,?);", -1, &estmt, NULL) != SQLITE_OK) {
        if (nstmt) sqlite3_finalize(nstmt);
        if (estmt) sqlite3_finalize(estmt);
        sql_exec(db, "ROLLBACK;");
        xfree(way_refs);
        idmap_free(&needed_nodes);
        fclose(fp);
        return -1;
    }

    size_t nodes_ins = 0, edges_ins = 0;
    in_way = 0;
    is_highway = 0;
    way_ref_cnt = 0;
    RoadClass rc = RC_UNCLASSIFIED;
    int oneway = 0; /* 0 = both, 1 = forward, -1 = reverse */
    float maxspd = -1.0f;
    float smooth = 0.85f;
    uint8_t wflags = EDGE_NONE;
    int motor_blocked = 0;

    LOG_INFO("OSM Import Pass 2: inserting highway nodes & edges into SQLite...");
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (!in_way) {
            if (strstr(line, "<node ") != NULL) {
                char id_s[64], lat_s[64], lon_s[64];
                if (xml_get_attr(line, "id", id_s, sizeof(id_s))) {
                    int64_t nid = strtoll(id_s, NULL, 10);
                    if (idmap_get(&needed_nodes, nid) != VMR_INVALID &&
                        xml_get_attr(line, "lat", lat_s, sizeof(lat_s)) &&
                        xml_get_attr(line, "lon", lon_s, sizeof(lon_s))) {
                        sqlite3_reset(nstmt);
                        sqlite3_bind_int64(nstmt, 1, nid);
                        sqlite3_bind_double(nstmt, 2, strtod(lat_s, NULL));
                        sqlite3_bind_double(nstmt, 3, strtod(lon_s, NULL));
                        sqlite3_step(nstmt);
                        nodes_ins++;
                    }
                }
            } else if (strstr(line, "<way ") != NULL) {
                in_way = 1;
                is_highway = 0;
                way_ref_cnt = 0;
                rc = RC_UNCLASSIFIED;
                oneway = 0;
                maxspd = -1.0f;
                smooth = 0.85f;
                wflags = EDGE_NONE;
                motor_blocked = 0;
                if (strstr(line, "/>") != NULL) in_way = 0;
            }
        } else {
            char key[128], val[128];
            if (strstr(line, "<nd ") != NULL && xml_get_attr(line, "ref", val, sizeof(val))) {
                if (way_ref_cnt < way_ref_cap) {
                    *(way_refs + way_ref_cnt++) = strtoll(val, NULL, 10);
                }
            } else if (strstr(line, "<tag ") != NULL &&
                       xml_get_attr(line, "k", key, sizeof(key)) &&
                       xml_get_attr(line, "v", val, sizeof(val))) {
                if (strcmp(key, "highway") == 0) {
                    if (parse_uk_highway_tag(val, &rc)) {
                        is_highway = 1;
                        if (rc == RC_MOTORWAY) oneway = 1;
                    }
                } else if (strcmp(key, "oneway") == 0) {
                    if (strcmp(val, "yes") == 0 || strcmp(val, "1") == 0 || strcmp(val, "true") == 0) {
                        oneway = 1;
                    } else if (strcmp(val, "-1") == 0 || strcmp(val, "reverse") == 0) {
                        oneway = -1;
                    } else if (strcmp(val, "no") == 0) {
                        oneway = 0;
                    }
                } else if (strcmp(key, "junction") == 0 && strcmp(val, "roundabout") == 0) {
                    oneway = 1;
                } else if (strcmp(key, "maxspeed") == 0) {
                    maxspd = parse_uk_maxspeed_kmh(val, rc);
                } else if (strcmp(key, "toll") == 0 && strcmp(val, "yes") == 0) {
                    wflags |= EDGE_TOLL;
                    smooth *= 0.6f;
                } else if (strcmp(key, "bridge") == 0 && strcmp(val, "yes") == 0) {
                    wflags |= EDGE_BRIDGE;
                } else if (strcmp(key, "surface") == 0) {
                    if (strcmp(val, "gravel") == 0 || strcmp(val, "dirt") == 0 || strcmp(val, "unpaved") == 0) {
                        smooth = 0.45f;
                    } else if (strcmp(val, "asphalt") == 0 || strcmp(val, "paved") == 0) {
                        smooth = 0.95f;
                    }
                } else if ((strcmp(key, "motor_vehicle") == 0 || strcmp(key, "access") == 0) &&
                           (strcmp(val, "no") == 0 || strcmp(val, "private") == 0)) {
                    motor_blocked = 1;
                }
            } else if (strstr(line, "</way>") != NULL) {
                if (is_highway && !motor_blocked && way_ref_cnt >= 2) {
                    float spd = (maxspd > 0.0f) ? maxspd : 0.0f;
                    if (spd <= 0.0f) {
                        float dummy = 0.0f;
                        road_class_defaults(rc, &spd, &dummy);
                    }
                    for (size_t i = 0; i + 1 < way_ref_cnt; i++) {
                        int64_t a = *(way_refs + i);
                        int64_t b = *(way_refs + i + 1);
                        if (oneway == 0 || oneway == 1) {
                            sqlite3_reset(estmt);
                            sqlite3_bind_int64(estmt, 1, a);
                            sqlite3_bind_int64(estmt, 2, b);
                            sqlite3_bind_int(estmt, 3, (int)rc);
                            sqlite3_bind_double(estmt, 4, (double)spd);
                            sqlite3_bind_double(estmt, 5, (double)smooth);
                            sqlite3_bind_int(estmt, 6, (int)wflags);
                            sqlite3_step(estmt);
                            edges_ins++;
                        }
                        if (oneway == 0 || oneway == -1) {
                            sqlite3_reset(estmt);
                            sqlite3_bind_int64(estmt, 1, b);
                            sqlite3_bind_int64(estmt, 2, a);
                            sqlite3_bind_int(estmt, 3, (int)rc);
                            sqlite3_bind_double(estmt, 4, (double)spd);
                            sqlite3_bind_double(estmt, 5, (double)smooth);
                            sqlite3_bind_int(estmt, 6, (int)wflags);
                            sqlite3_step(estmt);
                            edges_ins++;
                        }
                    }
                }
                in_way = 0;
            }
        }
    }

    sqlite3_finalize(nstmt);
    sqlite3_finalize(estmt);
    sql_exec(db, "COMMIT;");
    sql_exec(db, "PRAGMA synchronous = NORMAL;");

    xfree(way_refs);
    idmap_free(&needed_nodes);
    fclose(fp);
    LOG_INFO("OSM import finished: %zu nodes, %zu directed edges inserted", nodes_ins, edges_ins);
    return 0;
}

/* ========================================================================= */
/* SECTION 7B: OSM PBF (PROTOBUF) IMPORTER                                    */
/* ========================================================================= */

#define PBF_WAY_MAGIC 0x56525731u

typedef struct {
    const uint8_t *buf;
    size_t len;
    size_t pos;
    int err;
} PbReader;

static uint64_t pb_varint(PbReader *r) {
    uint64_t v = 0;
    int shift = 0;
    while (r->pos < r->len) {
        uint8_t b = *(r->buf + r->pos++);
        if (shift < 64) {
            v |= (uint64_t)(b & 0x7FU) << shift;
        }
        if ((b & 0x80U) == 0) {
            return v;
        }
        shift += 7;
        if (shift > 70) {
            r->err = 1;
            return 0;
        }
    }
    r->err = 1;
    return 0;
}

static int64_t pb_svarint(PbReader *r) {
    uint64_t v = pb_varint(r);
    return (int64_t)((v >> 1) ^ (0 - (v & 1U)));
}

static int pb_next(PbReader *r, int *field, int *wire) {
    if (r->pos >= r->len) {
        return 0;
    }
    uint64_t tag = pb_varint(r);
    if (r->err) {
        return 0;
    }
    *field = (int)(tag >> 3);
    *wire = (int)(tag & 7U);
    if (*field == 0) {
        r->err = 1;
        return 0;
    }
    return 1;
}

static const uint8_t *pb_bytes(PbReader *r, size_t *out_len) {
    uint64_t n = pb_varint(r);
    if (r->err || n > (uint64_t)(r->len - r->pos)) {
        r->err = 1;
        return NULL;
    }
    const uint8_t *p = r->buf + r->pos;
    r->pos += (size_t)n;
    *out_len = (size_t)n;
    return p;
}

static void pb_skip(PbReader *r, int wire) {
    switch (wire) {
        case 0: pb_varint(r); break;
        case 1:
            if (r->len - r->pos >= 8U) r->pos += 8U; else r->err = 1;
            break;
        case 2: { size_t n; (void)pb_bytes(r, &n); break; }
        case 5:
            if (r->len - r->pos >= 4U) r->pos += 4U; else r->err = 1;
            break;
        default: r->err = 1; break;
    }
}

typedef struct {
    const char *p;
    uint32_t n;
} PbStr;

typedef struct {
    PbStr *v;
    uint32_t n;
    uint32_t cap;
} PbStrTab;

typedef struct {
    const uint8_t *p;
    size_t n;
} PbRange;

typedef struct {
    const uint8_t *data;
    size_t size;
    PbStrTab strings;
    int64_t granularity;
    int64_t lat_off;
    int64_t lon_off;
} PbBlock;

typedef struct {
    int64_t  id;
    uint8_t  rc;
    uint8_t  flags;
    int8_t   oneway;
    float    spd;
    float    smooth;
    uint32_t ref_cnt;
} PbfWayRec;

typedef struct {
    int64_t  from_way;
    int64_t  to_way;
    int64_t  via_node;
    uint8_t  type;   /* 1=no_left 2=no_right 3=no_u_turn 4=no_straight */
} PbfTurnRes;

typedef struct {
    sqlite3 *db;
    sqlite3_stmt *nstmt;
    FILE *way_tmp;
    int pass;
    IdMap way_map;
    uint64_t *way_off;
    size_t way_off_n;
    size_t way_off_cap;
    PbfTurnRes *turn_res;
    size_t turn_res_n;
    size_t turn_res_cap;
    int64_t *needed;
    size_t needed_n;
    size_t needed_cap;
    int64_t *needed_sorted;
    size_t needed_sorted_n;
    size_t needed_idx;
    uint64_t ways_routable;
    uint64_t nodes_inserted;
    uint64_t data_blobs;
} PbfCtx;

static int pb_str_is(const PbStr *s, const char *lit) {
    if (!s || !s->p) {
        return 0;
    }
    size_t n = strlen(lit);
    if (s->n != (uint32_t)n) {
        return 0;
    }
    return n == 0 || memcmp(s->p, lit, n) == 0;
}

static int pb_str_cstr(const PbStr *s, char *out, size_t cap) {
    if (!s || !s->p || s->n >= cap) {
        return 0;
    }
    memcpy(out, s->p, s->n);
    *(out + s->n) = '\0';
    return 1;
}

static int pbf_strtab_add(PbStrTab *st, const uint8_t *p, size_t n) {
    if (st->n == st->cap) {
        uint32_t ncap = st->cap ? st->cap * 2U : 64U;
        PbStr *nv = (PbStr *)xrealloc(st->v, (size_t)ncap * sizeof(PbStr));
        if (!nv) {
            LOG_ERR("OOM: PBF string table");
            return -1;
        }
        st->v = nv;
        st->cap = ncap;
    }
    *(st->v + st->n) = (PbStr){ (const char *)p, (uint32_t)n };
    st->n++;
    return 0;
}

static void pbf_parse_string_table(PbStrTab *st, const uint8_t *p, size_t n) {
    PbReader r = { p, n, 0, 0 };
    int f, w;
    while (pb_next(&r, &f, &w)) {
        if (f == 1 && w == 2) {
            size_t l = 0;
            const uint8_t *s = pb_bytes(&r, &l);
            if (r.err) {
                return;
            }
            if (pbf_strtab_add(st, s, l) != 0) {
                return;
            }
        } else {
            pb_skip(&r, w);
        }
    }
}

static int cmp_i64(const void *a, const void *b) {
    int64_t x = *(const int64_t *)a;
    int64_t y = *(const int64_t *)b;
    return (x > y) - (x < y);
}

static int needed_contains(PbfCtx *c, int64_t id) {
    size_t n = c->needed_sorted_n;
    if (n == 0) {
        return 0;
    }
    size_t i = c->needed_idx;
    while (i < n && *(c->needed_sorted + i) < id) {
        i++;
    }
    c->needed_idx = i;
    if (i < n && *(c->needed_sorted + i) == id) {
        return 1;
    }
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (*(c->needed_sorted + mid) < id) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo < n && *(c->needed_sorted + lo) == id;
}

static void pbf_insert_node(PbfCtx *c, int64_t id, double lat, double lon) {
    sqlite3_stmt *st = c->nstmt;
    sqlite3_reset(st);
    sqlite3_bind_int64(st, 1, id);
    sqlite3_bind_double(st, 2, lat);
    sqlite3_bind_double(st, 3, lon);
    sqlite3_step(st);
    c->nodes_inserted++;
    if ((c->nodes_inserted % 1000000ULL) == 0) {
        LOG_INFO("PBF: %llu highway nodes inserted",
                 (unsigned long long)c->nodes_inserted);
    }
}

static void pbf_handle_node(PbfCtx *c, PbBlock *b, const uint8_t *p, size_t n) {
    PbReader r = { p, n, 0, 0 };
    int64_t id = 0, lat = 0, lon = 0;
    int have_id = 0, have_lat = 0, have_lon = 0;
    int f, w;
    while (pb_next(&r, &f, &w)) {
        if (f == 1 && w == 0) {
            id = (int64_t)pb_varint(&r);
            have_id = 1;
        } else if (f == 8 && w == 0) {
            lat = pb_svarint(&r);
            have_lat = 1;
        } else if (f == 9 && w == 0) {
            lon = pb_svarint(&r);
            have_lon = 1;
        } else {
            pb_skip(&r, w);
        }
    }
    if (r.err || !have_id || !have_lat || !have_lon || !needed_contains(c, id)) {
        return;
    }
    double la = ((double)b->lat_off + (double)b->granularity * (double)lat) / 1e9;
    double lo = ((double)b->lon_off + (double)b->granularity * (double)lon) / 1e9;
    pbf_insert_node(c, id, la, lo);
}

static void pbf_handle_dense(PbfCtx *c, PbBlock *b, const uint8_t *p, size_t n) {
    PbReader r = { p, n, 0, 0 };
    const uint8_t *idp = NULL, *latp = NULL, *lonp = NULL;
    size_t idn = 0, latn = 0, lonn = 0;
    int f, w;
    while (pb_next(&r, &f, &w)) {
        if (f == 1 && w == 2) {
            idp = pb_bytes(&r, &idn);
        } else if (f == 8 && w == 2) {
            latp = pb_bytes(&r, &latn);
        } else if (f == 9 && w == 2) {
            lonp = pb_bytes(&r, &lonn);
        } else {
            pb_skip(&r, w);
        }
    }
    if (r.err || !idp || !latp || !lonp) {
        return;
    }
    PbReader idr = { idp, idn, 0, 0 };
    PbReader lr = { latp, latn, 0, 0 };
    PbReader nr = { lonp, lonn, 0, 0 };
    int64_t id = 0, lat = 0, lon = 0;
    while (idr.pos < idr.len) {
        uint64_t a = pb_varint(&idr);
        uint64_t bb = pb_varint(&lr);
        uint64_t cc = pb_varint(&nr);
        if (idr.err || lr.err || nr.err) {
            break;
        }
        id += (int64_t)((a >> 1) ^ (0 - (a & 1U)));
        lat += (int64_t)((bb >> 1) ^ (0 - (bb & 1U)));
        lon += (int64_t)((cc >> 1) ^ (0 - (cc & 1U)));
        if (!needed_contains(c, id)) {
            continue;
        }
        double la = ((double)b->lat_off + (double)b->granularity * (double)lat) / 1e9;
        double lo = ((double)b->lon_off + (double)b->granularity * (double)lon) / 1e9;
        pbf_insert_node(c, id, la, lo);
    }
}

static void pbf_handle_way(PbfCtx *c, PbBlock *b, const uint8_t *p, size_t n) {
    PbReader r = { p, n, 0, 0 };
    int64_t id = 0;
    const uint8_t *keys = NULL, *vals = NULL, *refs = NULL;
    size_t kn = 0, vn = 0, rn = 0;
    int have_id = 0, have_keys = 0, have_vals = 0, have_refs = 0;
    int f, w;
    while (pb_next(&r, &f, &w)) {
        if (f == 1 && w == 0) {
            id = (int64_t)pb_varint(&r);
            have_id = 1;
        } else if (f == 2 && w == 2) {
            keys = pb_bytes(&r, &kn);
            have_keys = 1;
        } else if (f == 3 && w == 2) {
            vals = pb_bytes(&r, &vn);
            have_vals = 1;
        } else if (f == 8 && w == 2) {
            refs = pb_bytes(&r, &rn);
            have_refs = 1;
        } else {
            pb_skip(&r, w);
        }
    }
    if (r.err || !have_id || !have_refs) {
        return;
    }

    int64_t *ra = NULL;
    size_t rc = 0, rcap = 0;
    PbReader rr = { refs, rn, 0, 0 };
    int64_t cur = 0;
    while (rr.pos < rr.len) {
        uint64_t dv = pb_varint(&rr);
        if (rr.err) {
            break;
        }
        cur += (int64_t)((dv >> 1) ^ (0 - (dv & 1U)));
        if (rc == rcap) {
            size_t ncap = rcap ? rcap * 2U : 16U;
            int64_t *na = (int64_t *)xrealloc(ra, ncap * sizeof(int64_t));
            if (!na) {
                xfree(ra);
                LOG_ERR("OOM: PBF way refs");
                return;
            }
            ra = na;
            rcap = ncap;
        }
        *(ra + rc++) = cur;
    }
    if (rr.err || rc < 2) {
        xfree(ra);
        return;
    }

    RoadClass rclass = RC_UNCLASSIFIED;
    int is_highway = 0;
    int oneway = 0;
    int motor_blocked = 0;
    float maxspd = -1.0f;
    float smooth = 0.85f;
    uint8_t flags = EDGE_NONE;

    if (have_keys && have_vals) {
        PbReader kr = { keys, kn, 0, 0 };
        PbReader vr = { vals, vn, 0, 0 };
        while (kr.pos < kr.len && vr.pos < vr.len) {
            uint64_t ki = pb_varint(&kr);
            uint64_t vi = pb_varint(&vr);
            if (kr.err || vr.err) {
                break;
            }
            if (ki >= b->strings.n || vi >= b->strings.n) {
                continue;
            }
            const PbStr *k = b->strings.v + ki;
            const PbStr *v = b->strings.v + vi;
            char vb[128];
            if (pb_str_is(k, "highway")) {
                if (!pb_str_cstr(v, vb, sizeof vb)) {
                    continue;
                }
                RoadClass t;
                if (parse_uk_highway_tag(vb, &t)) {
                    rclass = t;
                    is_highway = 1;
                    if (rclass == RC_MOTORWAY) {
                        oneway = 1;
                    }
                }
            } else if (pb_str_is(k, "oneway")) {
                if (!pb_str_cstr(v, vb, sizeof vb)) {
                    continue;
                }
                if (strcmp(vb, "yes") == 0 || strcmp(vb, "1") == 0 || strcmp(vb, "true") == 0) {
                    oneway = 1;
                } else if (strcmp(vb, "-1") == 0 || strcmp(vb, "reverse") == 0) {
                    oneway = -1;
                } else if (strcmp(vb, "no") == 0) {
                    oneway = 0;
                }
            } else if (pb_str_is(k, "junction")) {
                if (pb_str_is(v, "roundabout")) {
                    oneway = 1;
                }
            } else if (pb_str_is(k, "maxspeed")) {
                if (!pb_str_cstr(v, vb, sizeof vb)) {
                    continue;
                }
                maxspd = parse_uk_maxspeed_kmh(vb, rclass);
            } else if (pb_str_is(k, "toll")) {
                if (pb_str_is(v, "yes")) {
                    flags |= EDGE_TOLL;
                    smooth *= 0.6f;
                }
            } else if (pb_str_is(k, "bridge")) {
                if (pb_str_is(v, "yes")) {
                    flags |= EDGE_BRIDGE;
                }
            } else if (pb_str_is(k, "surface")) {
                if (!pb_str_cstr(v, vb, sizeof vb)) {
                    continue;
                }
                if (strcmp(vb, "gravel") == 0 || strcmp(vb, "dirt") == 0 || strcmp(vb, "unpaved") == 0) {
                    smooth = 0.45f;
                } else if (strcmp(vb, "asphalt") == 0 || strcmp(vb, "paved") == 0) {
                    smooth = 0.95f;
                }
            } else if (pb_str_is(k, "motor_vehicle") || pb_str_is(k, "access")) {
                if (!pb_str_cstr(v, vb, sizeof vb)) {
                    continue;
                }
                if (strcmp(vb, "no") == 0 || strcmp(vb, "private") == 0) {
                    motor_blocked = 1;
                }
            }
        }
    }

    if (!is_highway || motor_blocked) {
        xfree(ra);
        return;
    }

    float spd = (maxspd > 0.0f) ? maxspd : 0.0f;
    if (spd <= 0.0f) {
        float dummy = 0.0f;
        road_class_defaults(rclass, &spd, &dummy);
    }

    for (size_t i = 0; i < rc; i++) {
        if (c->needed_n == c->needed_cap) {
            size_t ncap = c->needed_cap ? c->needed_cap * 2U : 1024U;
            int64_t *nn = (int64_t *)xrealloc(c->needed, ncap * sizeof(int64_t));
            if (!nn) {
                xfree(ra);
                LOG_ERR("OOM: PBF needed node set");
                return;
            }
            c->needed = nn;
            c->needed_cap = ncap;
        }
        *(c->needed + c->needed_n++) = *(ra + i);
    }

    PbfWayRec rec;
    memset(&rec, 0, sizeof rec);
    rec.id = id;
    rec.rc = (uint8_t)rclass;
    rec.flags = flags;
    rec.oneway = (int8_t)oneway;
    rec.spd = spd;
    rec.smooth = smooth;
    rec.ref_cnt = (uint32_t)rc;
    long off = ftell(c->way_tmp);
    if (off >= 0) {
        if (c->way_off_n == c->way_off_cap) {
            size_t ncap = c->way_off_cap ? c->way_off_cap * 2U : 1024U;
            uint64_t *no = (uint64_t *)xrealloc(c->way_off, ncap * sizeof(uint64_t));
            if (no) {
                c->way_off = no;
                c->way_off_cap = ncap;
            }
        }
        if (c->way_off_n < c->way_off_cap) {
            *(c->way_off + c->way_off_n) = (uint64_t)off;
            idmap_put(&c->way_map, rec.id, (uint32_t)c->way_off_n);
            c->way_off_n++;
        }
    }
    if (fwrite(&rec, sizeof rec, 1, c->way_tmp) != 1 ||
        (rc > 0 && fwrite(ra, sizeof(int64_t), rc, c->way_tmp) != rc)) {
        LOG_ERR("Failed writing PBF way cache");
    }
    c->ways_routable++;
    xfree(ra);
}

static void pbf_handle_relation(PbfCtx *c, PbBlock *b, const uint8_t *p, size_t n) {
    PbReader r = { p, n, 0, 0 };
    const uint8_t *keys = NULL, *vals = NULL, *roles = NULL, *memids = NULL, *types = NULL;
    size_t kn = 0, vn = 0, rolen = 0, mlen = 0, tlen = 0;
    int have_keys = 0, have_vals = 0, have_roles = 0, have_memids = 0, have_types = 0;
    int f, w;
    while (pb_next(&r, &f, &w)) {
        if (f == 2 && w == 2) { keys = pb_bytes(&r, &kn); have_keys = 1; }
        else if (f == 3 && w == 2) { vals = pb_bytes(&r, &vn); have_vals = 1; }
        else if (f == 8 && w == 2) { roles = pb_bytes(&r, &rolen); have_roles = 1; }
        else if (f == 9 && w == 2) { memids = pb_bytes(&r, &mlen); have_memids = 1; }
        else if (f == 10 && w == 2) { types = pb_bytes(&r, &tlen); have_types = 1; }
        else pb_skip(&r, w);
    }
    if (r.err || !have_roles || !have_memids || !have_types) {
        return;
    }

    int is_restriction = 0;
    int rtype = 0;
    int except_car = 0;
    if (have_keys && have_vals) {
        PbReader kr = { keys, kn, 0, 0 };
        PbReader vr = { vals, vn, 0, 0 };
        while (kr.pos < kr.len && vr.pos < vr.len) {
            uint64_t ki = pb_varint(&kr);
            uint64_t vi = pb_varint(&vr);
            if (kr.err || vr.err) break;
            if (ki >= b->strings.n || vi >= b->strings.n) continue;
            const PbStr *k = b->strings.v + ki;
            const PbStr *v = b->strings.v + vi;
            if (pb_str_is(k, "type")) {
                if (pb_str_is(v, "restriction")) is_restriction = 1;
            } else if (pb_str_is(k, "restriction")) {
                if (pb_str_is(v, "no_left_turn")) rtype = 1;
                else if (pb_str_is(v, "no_right_turn")) rtype = 2;
                else if (pb_str_is(v, "no_u_turn")) rtype = 3;
                else if (pb_str_is(v, "no_straight_on")) rtype = 4;
            } else if (pb_str_is(k, "except")) {
                if (pb_str_is(v, "motorcar") || pb_str_is(v, "car")) except_car = 1;
            }
        }
    }
    if (!is_restriction || rtype == 0 || except_car) {
        return;
    }

    PbReader rr = { roles, rolen, 0, 0 };
    PbReader mr = { memids, mlen, 0, 0 };
    PbReader tr = { types, tlen, 0, 0 };
    int64_t cur = 0;
    int64_t from_way = 0, to_way = 0, via_node = 0;
    int have_from = 0, have_to = 0, have_via = 0;
    while (rr.pos < rr.len) {
        uint64_t role = pb_varint(&rr);
        uint64_t dm = pb_varint(&mr);
        uint64_t ty = pb_varint(&tr);
        if (rr.err || mr.err || tr.err) break;
        cur += (int64_t)((dm >> 1) ^ (0 - (dm & 1U)));
        if (role >= b->strings.n) continue;
        const PbStr *rs = b->strings.v + role;
        if (ty == 1) {
            if (pb_str_is(rs, "from")) { from_way = cur; have_from = 1; }
            else if (pb_str_is(rs, "to")) { to_way = cur; have_to = 1; }
        } else if (ty == 0) {
            if (pb_str_is(rs, "via")) { via_node = cur; have_via = 1; }
        }
    }
    if (!have_from || !have_to || !have_via) {
        return;
    }
    if (c->turn_res_n == c->turn_res_cap) {
        size_t ncap = c->turn_res_cap ? c->turn_res_cap * 2U : 256U;
        PbfTurnRes *nr = (PbfTurnRes *)xrealloc(c->turn_res, ncap * sizeof(PbfTurnRes));
        if (!nr) {
            return;
        }
        c->turn_res = nr;
        c->turn_res_cap = ncap;
    }
    PbfTurnRes *out = c->turn_res + c->turn_res_n++;
    out->from_way = from_way;
    out->to_way = to_way;
    out->via_node = via_node;
    out->type = (uint8_t)rtype;
}

static int pbf_read_way_refs(PbfCtx *c, uint32_t idx, PbfWayRec *out_rec,
                             int64_t **out_refs, size_t *out_n) {
    if (idx >= c->way_off_n) {
        return -1;
    }
    if (fseek(c->way_tmp, (long)*(c->way_off + idx), SEEK_SET) != 0) {
        return -1;
    }
    PbfWayRec rec;
    if (fread(&rec, sizeof rec, 1, c->way_tmp) != 1) {
        return -1;
    }
    int64_t *refs = NULL;
    if (rec.ref_cnt > 0) {
        refs = (int64_t *)xmalloc((size_t)rec.ref_cnt * sizeof(int64_t));
        if (!refs || fread(refs, sizeof(int64_t), rec.ref_cnt, c->way_tmp) != rec.ref_cnt) {
            xfree(refs);
            return -1;
        }
    }
    if (rec.oneway == -1) {
        for (uint32_t i = 0; i < rec.ref_cnt / 2U; i++) {
            int64_t tmp = *(refs + i);
            *(refs + i) = *(refs + (rec.ref_cnt - 1U - i));
            *(refs + (rec.ref_cnt - 1U - i)) = tmp;
        }
    }
    *out_rec = rec;
    *out_refs = refs;
    *out_n = rec.ref_cnt;
    return 0;
}

static int pbf_apply_restrictions(PbfCtx *c) {
    if (c->turn_res_n == 0) {
        LOG_INFO("PBF: no turn restrictions found");
        return 0;
    }
    sqlite3_stmt *st = NULL;
    const char *sql = "INSERT OR REPLACE INTO turn_restrictions(via_node,from_node,to_node,mode) VALUES(?,?,?,?);";
    if (sqlite3_prepare_v2(c->db, sql, -1, &st, NULL) != SQLITE_OK) {
        LOG_ERR("Turn restriction prepare: %s", sqlite3_errmsg(c->db));
        return -1;
    }
    if (sql_exec(c->db, "BEGIN") != 0) {
        sqlite3_finalize(st);
        return -1;
    }
    uint64_t applied = 0;
    uint64_t skipped = 0;
    for (size_t i = 0; i < c->turn_res_n; i++) {
        PbfTurnRes *rr = c->turn_res + i;
        uint32_t fi = idmap_get(&c->way_map, rr->from_way);
        uint32_t ti = idmap_get(&c->way_map, rr->to_way);
        if (fi == VMR_INVALID || ti == VMR_INVALID) {
            skipped++;
            continue;
        }
        PbfWayRec fr, tr;
        int64_t *frefs = NULL, *trefs = NULL;
        size_t fn = 0, tn = 0;
        if (pbf_read_way_refs(c, fi, &fr, &frefs, &fn) != 0 ||
            pbf_read_way_refs(c, ti, &tr, &trefs, &tn) != 0) {
            xfree(frefs);
            xfree(trefs);
            skipped++;
            continue;
        }
        size_t kf = SIZE_MAX, kt = SIZE_MAX;
        for (size_t k = 0; k < fn; k++) {
            if (*(frefs + k) == rr->via_node) { kf = k; break; }
        }
        for (size_t k = 0; k < tn; k++) {
            if (*(trefs + k) == rr->via_node) { kt = k; break; }
        }
        if (kf == SIZE_MAX || kf == 0 || kt == SIZE_MAX || kt + 1 >= tn) {
            xfree(frefs);
            xfree(trefs);
            skipped++;
            continue;
        }
        sqlite3_reset(st);
        sqlite3_bind_int64(st, 1, rr->via_node);
        sqlite3_bind_int64(st, 2, *(frefs + kf - 1));
        sqlite3_bind_int64(st, 3, *(trefs + kt + 1));
        sqlite3_bind_int(st, 4, 0);
        sqlite3_step(st);
        applied++;
        xfree(frefs);
        xfree(trefs);
    }
    sql_exec(c->db, "COMMIT");
    sqlite3_finalize(st);
    LOG_INFO("PBF: applied %llu turn restrictions (%llu skipped)",
             (unsigned long long)applied, (unsigned long long)skipped);
    return 0;
}

static int pbf_parse_group(PbfCtx *c, PbBlock *b, const uint8_t *p, size_t n) {
    PbReader r = { p, n, 0, 0 };
    int f, w;
    while (pb_next(&r, &f, &w)) {
        if (w != 2) {
            pb_skip(&r, w);
            continue;
        }
        size_t l = 0;
        const uint8_t *q = pb_bytes(&r, &l);
        if (r.err) {
            return -1;
        }
        if (f == 1) {
            if (c->pass == 2) {
                pbf_handle_node(c, b, q, l);
            }
        } else if (f == 2) {
            if (c->pass == 2) {
                pbf_handle_dense(c, b, q, l);
            }
        } else if (f == 3) {
            if (c->pass == 1) {
                pbf_handle_way(c, b, q, l);
            }
        } else if (f == 4) {
            if (c->pass == 1) {
                pbf_handle_relation(c, b, q, l);
            }
        }
    }
    return r.err ? -1 : 0;
}

static int pbf_parse_primitive_block(PbfCtx *c, const uint8_t *data, size_t size) {
    PbBlock b;
    memset(&b, 0, sizeof b);
    b.data = data;
    b.size = size;
    b.granularity = 100;
    b.lat_off = 0;
    b.lon_off = 0;

    PbRange *groups = NULL;
    int ng = 0, gcap = 0;
    PbReader r = { data, size, 0, 0 };
    int f, w;
    while (pb_next(&r, &f, &w)) {
        if (f == 1 && w == 2) {
            size_t l = 0;
            const uint8_t *s = pb_bytes(&r, &l);
            if (r.err) {
                break;
            }
            pbf_parse_string_table(&b.strings, s, l);
        } else if (f == 2 && w == 2) {
            size_t l = 0;
            const uint8_t *s = pb_bytes(&r, &l);
            if (r.err) {
                break;
            }
            if (ng == gcap) {
                int ncap = gcap ? gcap * 2 : 8;
                PbRange *tmp = (PbRange *)xrealloc(groups, (size_t)ncap * sizeof(PbRange));
                if (!tmp) {
                    LOG_ERR("OOM: PBF primitive groups");
                    xfree(groups);
                    xfree(b.strings.v);
                    return -1;
                }
                groups = tmp;
                gcap = ncap;
            }
            *(groups + ng) = (PbRange){ s, l };
            ng++;
        } else if (f == 3 && w == 0) {
            b.granularity = (int64_t)pb_varint(&r);
        } else if (f == 5 && w == 0) {
            b.lat_off = (int64_t)pb_varint(&r);
        } else if (f == 6 && w == 0) {
            b.lon_off = (int64_t)pb_varint(&r);
        } else {
            pb_skip(&r, w);
        }
    }

    int rc = r.err ? -1 : 0;
    if (rc == 0) {
        for (int i = 0; i < ng; i++) {
            if (pbf_parse_group(c, &b, groups[i].p, groups[i].n) != 0) {
                rc = -1;
                break;
            }
        }
    }
    xfree(groups);
    xfree(b.strings.v);
    return rc;
}

static int pbf_process_blob(PbfCtx *c, const uint8_t *blob, size_t blob_len) {
    PbReader r = { blob, blob_len, 0, 0 };
    const uint8_t *raw = NULL, *zdata = NULL;
    size_t raw_len = 0, z_len = 0;
    uint64_t raw_size = 0;
    int have_raw_size = 0;
    int f, w;
    while (pb_next(&r, &f, &w)) {
        if (f == 1 && w == 2) {
            raw = pb_bytes(&r, &raw_len);
        } else if (f == 2 && w == 0) {
            raw_size = pb_varint(&r);
            have_raw_size = 1;
        } else if (f == 3 && w == 2) {
            zdata = pb_bytes(&r, &z_len);
        } else if (f == 4 && w == 2) {
            LOG_ERR("LZMA-compressed PBF blobs are not supported");
            return -1;
        } else if (f == 7 && w == 2) {
            LOG_ERR("ZSTD-compressed PBF blobs are not supported");
            return -1;
        } else {
            pb_skip(&r, w);
        }
    }
    if (r.err) {
        return -1;
    }

    const uint8_t *payload = raw;
    size_t payload_len = raw_len;
    uint8_t *inflated = NULL;
    if (zdata) {
        if (!have_raw_size || raw_size == 0 || raw_size > (128ULL << 20)) {
            LOG_ERR("Invalid PBF raw_size");
            return -1;
        }
        inflated = (uint8_t *)xmalloc((size_t)raw_size);
        if (!inflated) {
            return -1;
        }
        uLongf dl = (uLongf)raw_size;
        int zr = uncompress(inflated, &dl, zdata, (uLong)z_len);
        if (zr != Z_OK) {
            LOG_ERR("zlib inflate failed (%d)", zr);
            xfree(inflated);
            return -1;
        }
        payload = inflated;
        payload_len = (size_t)dl;
    } else if (!raw) {
        LOG_ERR("PBF blob has no payload");
        return -1;
    }

    int rc = pbf_parse_primitive_block(c, payload, payload_len);
    xfree(inflated);
    return rc;
}

static int pbf_read_blob_header(FILE *fp, char *type, size_t type_cap, uint32_t *datasize) {
    uint8_t lb[4];
    size_t got = fread(lb, 1, 4, fp);
    if (got == 0) {
        return 0;
    }
    if (got != 4) {
        LOG_ERR("Truncated PBF blob header length");
        return -1;
    }
    uint32_t len = ((uint32_t)lb[0] << 24) | ((uint32_t)lb[1] << 16) |
                   ((uint32_t)lb[2] << 8) | (uint32_t)lb[3];
    if (len == 0 || len > (1U << 20)) {
        LOG_ERR("Invalid PBF blob header length %u", len);
        return -1;
    }
    uint8_t *hdr = (uint8_t *)xmalloc(len);
    if (!hdr) {
        return -1;
    }
    if (fread(hdr, 1, len, fp) != len) {
        xfree(hdr);
        LOG_ERR("Truncated PBF blob header");
        return -1;
    }
    PbReader r = { hdr, len, 0, 0 };
    type[0] = '\0';
    *datasize = 0;
    int f, w;
    while (pb_next(&r, &f, &w)) {
        if (f == 1 && w == 2) {
            size_t l = 0;
            const uint8_t *s = pb_bytes(&r, &l);
            if (r.err) {
                break;
            }
            size_t n = (l < type_cap - 1U) ? l : (type_cap - 1U);
            memcpy(type, s, n);
            *(type + n) = '\0';
        } else if (f == 3 && w == 0) {
            *datasize = (uint32_t)pb_varint(&r);
        } else {
            pb_skip(&r, w);
        }
    }
    xfree(hdr);
    if (r.err || *datasize == 0) {
        LOG_ERR("Invalid PBF blob header");
        return -1;
    }
    return 1;
}

static int pbf_skip_bytes(FILE *fp, uint32_t n) {
    uint8_t buf[65536];
    while (n > 0) {
        size_t chunk = (n < (uint32_t)sizeof(buf)) ? n : (size_t)sizeof(buf);
        if (fread(buf, 1, chunk, fp) != chunk) {
            return -1;
        }
        n -= (uint32_t)chunk;
    }
    return 0;
}

static int pbf_parse_file(PbfCtx *c, FILE *fp) {
    uint8_t *blob = NULL;
    size_t cap = 0;
    uint64_t blobs = 0;
    for (;;) {
        char type[64];
        uint32_t dsize = 0;
        int rr = pbf_read_blob_header(fp, type, sizeof type, &dsize);
        if (rr == 0) {
            break;
        }
        if (rr < 0) {
            xfree(blob);
            return -1;
        }
        blobs++;
        if (strcmp(type, "OSMData") != 0) {
            if (pbf_skip_bytes(fp, dsize) != 0) {
                xfree(blob);
                LOG_ERR("Truncated PBF file");
                return -1;
            }
            continue;
        }
        if (cap < dsize) {
            uint8_t *nb = (uint8_t *)xrealloc(blob, dsize);
            if (!nb) {
                xfree(blob);
                return -1;
            }
            blob = nb;
            cap = dsize;
        }
        if (fread(blob, 1, dsize, fp) != dsize) {
            xfree(blob);
            LOG_ERR("Truncated PBF data blob");
            return -1;
        }
        c->data_blobs++;
        if (pbf_process_blob(c, blob, dsize) != 0) {
            xfree(blob);
            return -1;
        }
        if ((c->data_blobs & 0xFFU) == 0) {
            LOG_INFO("PBF pass %d: %llu data blobs", c->pass,
                     (unsigned long long)c->data_blobs);
        }
    }
    xfree(blob);
    LOG_INFO("PBF pass %d complete: %llu blobs", c->pass, (unsigned long long)blobs);
    return 0;
}

static int pbf_insert_edges(sqlite3 *db, FILE *fp) {
    uint32_t magic = 0;
    if (fread(&magic, sizeof magic, 1, fp) != 1 || magic != PBF_WAY_MAGIC) {
        LOG_ERR("Invalid PBF way cache");
        return -1;
    }
    sqlite3_stmt *st = NULL;
    const char *sql =
        "INSERT INTO osm_edges(from_node,to_node,road_class,speed_kmh,smooth,flags) "
        "VALUES(?,?,?,?,?,?);";
    if (sqlite3_prepare_v2(db, sql, -1, &st, NULL) != SQLITE_OK) {
        LOG_ERR("Edge insert prepare: %s", sqlite3_errmsg(db));
        return -1;
    }
    if (sql_exec(db, "BEGIN") != 0) {
        sqlite3_finalize(st);
        return -1;
    }

    PbfWayRec rec;
    uint64_t ways = 0, edges = 0;
    while (fread(&rec, sizeof rec, 1, fp) == 1) {
        int64_t *refs = NULL;
        if (rec.ref_cnt > 0) {
            refs = (int64_t *)xmalloc((size_t)rec.ref_cnt * sizeof(int64_t));
            if (!refs || fread(refs, sizeof(int64_t), rec.ref_cnt, fp) != rec.ref_cnt) {
                xfree(refs);
                LOG_ERR("Truncated PBF way cache");
                break;
            }
        }
        for (uint32_t i = 0; i + 1 < rec.ref_cnt; i++) {
            int64_t a = *(refs + i);
            int64_t bb = *(refs + i + 1);
            if (rec.oneway == 0 || rec.oneway == 1) {
                sqlite3_reset(st);
                sqlite3_bind_int64(st, 1, a);
                sqlite3_bind_int64(st, 2, bb);
                sqlite3_bind_int(st, 3, (int)rec.rc);
                sqlite3_bind_double(st, 4, (double)rec.spd);
                sqlite3_bind_double(st, 5, (double)rec.smooth);
                sqlite3_bind_int(st, 6, (int)rec.flags);
                sqlite3_step(st);
                edges++;
            }
            if (rec.oneway == 0 || rec.oneway == -1) {
                sqlite3_reset(st);
                sqlite3_bind_int64(st, 1, bb);
                sqlite3_bind_int64(st, 2, a);
                sqlite3_bind_int(st, 3, (int)rec.rc);
                sqlite3_bind_double(st, 4, (double)rec.spd);
                sqlite3_bind_double(st, 5, (double)rec.smooth);
                sqlite3_bind_int(st, 6, (int)rec.flags);
                sqlite3_step(st);
                edges++;
            }
        }
        xfree(refs);
        ways++;
        if ((ways & 0x3FFFFU) == 0) {
            LOG_INFO("PBF: %llu ways, %llu directed edges",
                     (unsigned long long)ways, (unsigned long long)edges);
        }
    }
    sql_exec(db, "COMMIT");
    sqlite3_finalize(st);
    LOG_INFO("PBF edge pass complete: %llu ways, %llu directed edges",
             (unsigned long long)ways, (unsigned long long)edges);
    return 0;
}

static int db_import_osm_pbf(sqlite3 *db, const char *pbf_path, int restrictions_only) {
    FILE *fp = fopen(pbf_path, "rb");
    if (!fp) {
        LOG_ERR("Cannot open PBF '%s'", pbf_path);
        return -1;
    }

    PbfCtx c;
    memset(&c, 0, sizeof c);
    c.db = db;
    c.pass = 1;
    c.way_tmp = tmpfile();
    char *tmp_path = NULL;
    if (!c.way_tmp) {
        size_t n = strlen(pbf_path) + 16U;
        tmp_path = (char *)xmalloc(n);
        if (tmp_path) {
            snprintf(tmp_path, n, "%s.ways.tmp", pbf_path);
            c.way_tmp = fopen(tmp_path, "wb+");
        }
        if (!c.way_tmp) {
            LOG_ERR("Cannot create temporary PBF way cache");
            xfree(tmp_path);
            fclose(fp);
            return -1;
        }
    }

    uint32_t magic = PBF_WAY_MAGIC;
    fwrite(&magic, sizeof magic, 1, c.way_tmp);

    if (sqlite3_prepare_v2(db,
            "INSERT OR IGNORE INTO osm_nodes(id,lat,lon) VALUES(?,?,?);",
            -1, &c.nstmt, NULL) != SQLITE_OK) {
        LOG_ERR("Node insert prepare: %s", sqlite3_errmsg(db));
        sqlite3_finalize(c.nstmt);
        fclose(c.way_tmp);
        fclose(fp);
        if (tmp_path) {
            xfree(tmp_path);
        }
        return -1;
    }

    sql_exec(db, "PRAGMA journal_mode = OFF;");
    sql_exec(db, "PRAGMA synchronous = OFF;");
    sql_exec(db, "PRAGMA locking_mode = EXCLUSIVE;");

    LOG_INFO("PBF pass 1/2: collecting routable ways from '%s'", pbf_path);
    int rc = pbf_parse_file(&c, fp);

    if (rc == 0 && c.needed_n > 0) {
        qsort(c.needed, c.needed_n, sizeof(int64_t), cmp_i64);
        size_t u = 1;
        for (size_t i = 1; i < c.needed_n; i++) {
            if (*(c.needed + i) != *(c.needed + u - 1U)) {
                *(c.needed + u++) = *(c.needed + i);
            }
        }
        c.needed_n = u;
    }
    c.needed_sorted = c.needed;
    c.needed_sorted_n = c.needed_n;
    c.needed_idx = 0;
    LOG_INFO("PBF pass 1 complete: %llu routable ways, %llu unique highway nodes",
             (unsigned long long)c.ways_routable, (unsigned long long)c.needed_sorted_n);

    if (rc == 0) {
        if (pbf_apply_restrictions(&c) != 0) {
            rc = -1;
        }
    }

    if (rc == 0 && !restrictions_only && c.needed_sorted_n > 0) {
        c.pass = 2;
        clearerr(fp);
        if (fseek(fp, 0, SEEK_SET) != 0) {
            rc = -1;
        }
        if (rc == 0 && sql_exec(db, "BEGIN") != 0) {
            rc = -1;
        }
        if (rc == 0) {
            LOG_INFO("PBF pass 2/2: inserting highway nodes");
            rc = pbf_parse_file(&c, fp);
        }
        if (rc == 0) {
            sql_exec(db, "COMMIT");
        } else {
            sql_exec(db, "ROLLBACK");
        }
    }

    if (rc == 0 && !restrictions_only) {
        fflush(c.way_tmp);
        rewind(c.way_tmp);
        if (pbf_insert_edges(db, c.way_tmp) != 0) {
            rc = -1;
        }
    }

    if (c.nstmt) {
        sqlite3_finalize(c.nstmt);
    }
    fclose(c.way_tmp);
    fclose(fp);
    if (tmp_path) {
        remove(tmp_path);
        xfree(tmp_path);
    }
    idmap_free(&c.way_map);
    xfree(c.way_off);
    xfree(c.turn_res);
    xfree(c.needed);

    if (rc != 0) {
        LOG_ERR("PBF import failed");
        return -1;
    }
    if (restrictions_only) {
        LOG_INFO("PBF restrictions-only pass complete: %llu ways, %llu data blobs",
                 (unsigned long long)c.ways_routable, (unsigned long long)c.data_blobs);
    } else {
        LOG_INFO("PBF import complete: %llu nodes, %llu ways, %llu data blobs",
                 (unsigned long long)c.nodes_inserted,
                 (unsigned long long)c.ways_routable,
                 (unsigned long long)c.data_blobs);
    }
    return 0;
}

static int db_import_csv(sqlite3 *db, const char *nodes_csv, const char *edges_csv) {
    FILE *fn = fopen(nodes_csv, "r");
    FILE *fe = fopen(edges_csv, "r");
    if (!fn || !fe) {
        if (fn) fclose(fn);
        if (fe) fclose(fe);
        LOG_ERR("Unable to open CSV inputs '%s' / '%s'", nodes_csv, edges_csv);
        return -1;
    }

    sql_exec(db, "PRAGMA synchronous = OFF;");
    sql_exec(db, "BEGIN TRANSACTION;");

    sqlite3_stmt *nstmt = NULL, *estmt = NULL;
    sqlite3_prepare_v2(db, "INSERT OR IGNORE INTO osm_nodes(id,lat,lon) VALUES(?,?,?);", -1, &nstmt, NULL);
    sqlite3_prepare_v2(db, "INSERT INTO osm_edges(from_node,to_node,road_class,speed_kmh,smooth,flags) VALUES(?,?,?,?,?,?);",
                       -1, &estmt, NULL);

    char line[512];
    size_t nc = 0, ec = 0;
    while (fgets(line, sizeof(line), fn) != NULL) {
        int64_t id = 0;
        double lat = 0.0, lon = 0.0;
        if (sscanf(line, "%lld,%lf,%lf", (long long *)&id, &lat, &lon) == 3) {
            sqlite3_reset(nstmt);
            sqlite3_bind_int64(nstmt, 1, id);
            sqlite3_bind_double(nstmt, 2, lat);
            sqlite3_bind_double(nstmt, 3, lon);
            sqlite3_step(nstmt);
            nc++;
        }
    }
    while (fgets(line, sizeof(line), fe) != NULL) {
        int64_t u = 0, v = 0;
        int rc = 4, flags = 0;
        double spd = 50.0, sm = 0.85;
        int parsed = sscanf(line, "%lld,%lld,%d,%lf,%lf,%d",
                            (long long *)&u, (long long *)&v, &rc, &spd, &sm, &flags);
        if (parsed >= 2) {
            sqlite3_reset(estmt);
            sqlite3_bind_int64(estmt, 1, u);
            sqlite3_bind_int64(estmt, 2, v);
            sqlite3_bind_int(estmt, 3, rc);
            sqlite3_bind_double(estmt, 4, spd);
            sqlite3_bind_double(estmt, 5, sm);
            sqlite3_bind_int(estmt, 6, flags);
            sqlite3_step(estmt);
            ec++;
        }
    }
    sqlite3_finalize(nstmt);
    sqlite3_finalize(estmt);
    sql_exec(db, "COMMIT;");
    sql_exec(db, "PRAGMA synchronous = NORMAL;");
    fclose(fn);
    fclose(fe);
    LOG_INFO("CSV import finished: %zu nodes, %zu edges", nc, ec);
    return 0;
}

/*
 * Seeds a realistic UK-wide inter-city test network into SQLite (London, Dover,
 * Southampton, Bristol, Cardiff, Birmingham, Norwich, Sheffield, Manchester,
 * Liverpool, Leeds, Hull, Newcastle, Edinburgh, Glasgow, Aberdeen, Inverness,
 * plus curved intermediate shape points to exercise degree-2 contraction,
 * cul-de-sac spurs to exercise 1-core peeling, and Isle of Wight / Belfast
 * to exercise disconnected component checking).
 */
typedef struct {
    int64_t id;
    double  lat;
    double  lon;
} DemoCity;

static int db_seed_uk_demo(sqlite3 *db) {
    sqlite3_stmt *chk = NULL;
    int existing = 0;
    if (sqlite3_prepare_v2(db, "SELECT COUNT(*) FROM osm_nodes;", -1, &chk, NULL) == SQLITE_OK) {
        if (sqlite3_step(chk) == SQLITE_ROW) existing = sqlite3_column_int(chk, 0);
        sqlite3_finalize(chk);
    }
    if (existing > 0) {
        LOG_INFO("Database already contains %d nodes; skipping UK demo seed.", existing);
        return 0;
    }

    static const DemoCity cities[] = {
        {  1, 51.5074, -0.1278 }, /* London */
        {  2, 51.1279,  1.3134 }, /* Dover */
        {  3, 50.9097, -1.4044 }, /* Southampton */
        {  4, 50.3755, -4.1427 }, /* Plymouth */
        {  5, 51.4545, -2.5879 }, /* Bristol */
        {  6, 51.4816, -3.1791 }, /* Cardiff */
        {  7, 51.7520, -1.2577 }, /* Oxford */
        {  8, 52.6309,  1.2974 }, /* Norwich */
        {  9, 52.2053,  0.1218 }, /* Cambridge */
        { 10, 52.4862, -1.8904 }, /* Birmingham */
        { 11, 52.9548, -1.1581 }, /* Nottingham */
        { 12, 53.3811, -1.4701 }, /* Sheffield */
        { 13, 53.4808, -2.2426 }, /* Manchester */
        { 14, 53.4084, -2.9916 }, /* Liverpool */
        { 15, 53.8008, -1.5491 }, /* Leeds */
        { 16, 53.7676, -0.3274 }, /* Hull */
        { 17, 54.8925, -2.9329 }, /* Carlisle */
        { 18, 54.9783, -1.6178 }, /* Newcastle */
        { 19, 55.9533, -3.1883 }, /* Edinburgh */
        { 20, 55.8642, -4.2518 }, /* Glasgow */
        { 21, 57.1497, -2.0943 }, /* Aberdeen */
        { 22, 57.4778, -4.2247 }, /* Inverness */
        /* Disconnected islands for island component test */
        { 30, 50.6938, -1.3047 }, /* Isle of Wight (Newport) */
        { 31, 50.7610, -1.2980 }, /* Isle of Wight (Cowes) */
        { 32, 54.5973, -5.9301 }, /* Belfast */
        { 33, 54.9966, -7.3086 }  /* Derry */
    };
    const size_t ncities = sizeof(cities) / sizeof(cities[0]);

    typedef struct {
        int64_t a, b;
        RoadClass rc;
        double spd;
        double smooth;
        int flags;
        int way_pts;
    } Corridor;

    static const Corridor corridors[] = {
        {  1,  2, RC_MOTORWAY,  112.6, 0.92, EDGE_NONE,   6 }, /* M20 London-Dover */
        {  1,  3, RC_MOTORWAY,  112.6, 0.90, EDGE_NONE,   5 }, /* M3 London-Southampton */
        {  1,  7, RC_MOTORWAY,  112.6, 0.92, EDGE_NONE,   4 }, /* M40 London-Oxford */
        {  1,  9, RC_MOTORWAY,  112.6, 0.90, EDGE_NONE,   4 }, /* M11 London-Cambridge */
        {  1, 10, RC_MOTORWAY,  112.6, 0.93, EDGE_NONE,   7 }, /* M1/M6 London-Birmingham */
        {  1,  5, RC_MOTORWAY,  112.6, 0.92, EDGE_NONE,   6 }, /* M4 London-Bristol */
        {  3,  4, RC_TRUNK,      96.5, 0.84, EDGE_NONE,   8 }, /* A35/A38 Soton-Plymouth */
        {  3,  5, RC_PRIMARY,    80.5, 0.80, EDGE_NONE,   5 }, /* A36 Soton-Bristol */
        {  4,  5, RC_MOTORWAY,  112.6, 0.91, EDGE_NONE,   7 }, /* M5/A38 Plymouth-Bristol */
        {  5,  6, RC_MOTORWAY,  112.6, 0.88, EDGE_BRIDGE, 3 }, /* M4 Severn Bridge */
        {  5,  7, RC_PRIMARY,    80.5, 0.82, EDGE_NONE,   4 }, /* A420 Bristol-Oxford */
        {  5, 10, RC_MOTORWAY,  112.6, 0.91, EDGE_NONE,   5 }, /* M5 Bristol-Birmingham */
        {  6, 10, RC_TRUNK,      96.5, 0.83, EDGE_NONE,   6 }, /* A449/M50 Cardiff-Birmingham */
        {  7, 10, RC_MOTORWAY,  112.6, 0.93, EDGE_NONE,   4 }, /* M40 Oxford-Birmingham */
        {  7,  9, RC_SECONDARY,  64.4, 0.78, EDGE_NONE,   5 }, /* A421 Oxford-Cambridge */
        {  9,  8, RC_TRUNK,      96.5, 0.86, EDGE_NONE,   5 }, /* A11 Cambridge-Norwich */
        {  9, 11, RC_TRUNK,      96.5, 0.85, EDGE_NONE,   5 }, /* A14/A1 Cambridge-Notts */
        {  8, 16, RC_PRIMARY,    80.5, 0.78, EDGE_NONE,   8 }, /* A17/A15 Norwich-Hull */
        { 10, 11, RC_MOTORWAY,  112.6, 0.91, EDGE_NONE,   4 }, /* M42/M1 Bham-Notts */
        { 10, 13, RC_MOTORWAY,  112.6, 0.92, EDGE_NONE,   6 }, /* M6 Bham-Manchester */
        { 10, 14, RC_MOTORWAY,  112.6, 0.89, EDGE_TOLL,   6 }, /* M6 Toll Bham-Liverpool */
        { 11, 12, RC_MOTORWAY,  112.6, 0.91, EDGE_NONE,   4 }, /* M1 Notts-Sheffield */
        { 12, 13, RC_PRIMARY,    72.4, 0.72, EDGE_NONE,   6 }, /* A57 Snake Pass */
        { 12, 15, RC_MOTORWAY,  112.6, 0.92, EDGE_NONE,   3 }, /* M1 Sheffield-Leeds */
        { 12, 16, RC_MOTORWAY,  112.6, 0.90, EDGE_BRIDGE, 4 }, /* M18/M62/Humber */
        { 13, 14, RC_MOTORWAY,  112.6, 0.94, EDGE_NONE,   3 }, /* M62 Manchester-Liverpool */
        { 13, 15, RC_MOTORWAY,  112.6, 0.90, EDGE_NONE,   4 }, /* M62 Manchester-Leeds */
        { 13, 17, RC_MOTORWAY,  112.6, 0.93, EDGE_NONE,   8 }, /* M6 Manchester-Carlisle */
        { 15, 16, RC_MOTORWAY,  112.6, 0.91, EDGE_NONE,   4 }, /* M62 Leeds-Hull */
        { 15, 18, RC_MOTORWAY,  112.6, 0.93, EDGE_NONE,   6 }, /* A1(M) Leeds-Newcastle */
        { 16, 18, RC_PRIMARY,    80.5, 0.80, EDGE_NONE,   7 }, /* A165/A19 Hull-Newcastle */
        { 17, 18, RC_TRUNK,      96.5, 0.85, EDGE_NONE,   4 }, /* A69 Carlisle-Newcastle */
        { 17, 20, RC_MOTORWAY,  112.6, 0.94, EDGE_NONE,   6 }, /* M74 Carlisle-Glasgow */
        { 17, 19, RC_PRIMARY,    80.5, 0.79, EDGE_NONE,   6 }, /* A7 Carlisle-Edinburgh */
        { 18, 19, RC_TRUNK,      96.5, 0.86, EDGE_NONE,   7 }, /* A1 Newcastle-Edinburgh */
        { 19, 20, RC_MOTORWAY,  112.6, 0.94, EDGE_NONE,   4 }, /* M8 Edinburgh-Glasgow */
        { 19, 21, RC_TRUNK,      96.5, 0.86, EDGE_BRIDGE, 8 }, /* M90/A90 Forth Bridge-Aberdeen */
        { 19, 22, RC_TRUNK,      96.5, 0.84, EDGE_NONE,   9 }, /* M9/A9 Edinburgh-Inverness */
        { 20, 22, RC_PRIMARY,    80.5, 0.77, EDGE_NONE,  10 }, /* A82 Glasgow-Fort William-Inverness */
        { 21, 22, RC_TRUNK,      88.5, 0.82, EDGE_NONE,   7 }, /* A96 Aberdeen-Inverness */
        /* Island links */
        { 30, 31, RC_PRIMARY,    64.4, 0.85, EDGE_NONE,   2 }, /* Isle of Wight */
        { 32, 33, RC_TRUNK,      96.5, 0.88, EDGE_NONE,   4 }  /* M2/A6 Belfast-Derry */
    };
    const size_t ncor = sizeof(corridors) / sizeof(corridors[0]);

    sql_exec(db, "BEGIN TRANSACTION;");
    sqlite3_stmt *nstmt = NULL, *estmt = NULL;
    sqlite3_prepare_v2(db, "INSERT INTO osm_nodes(id,lat,lon) VALUES(?,?,?);", -1, &nstmt, NULL);
    sqlite3_prepare_v2(db, "INSERT INTO osm_edges(from_node,to_node,road_class,speed_kmh,smooth,flags) VALUES(?,?,?,?,?,?);",
                       -1, &estmt, NULL);

    for (size_t i = 0; i < ncities; i++) {
        sqlite3_reset(nstmt);
        sqlite3_bind_int64(nstmt, 1, cities[i].id);
        sqlite3_bind_double(nstmt, 2, cities[i].lat);
        sqlite3_bind_double(nstmt, 3, cities[i].lon);
        sqlite3_step(nstmt);
    }

    int64_t next_id = 1000;
    for (size_t c = 0; c < ncor; c++) {
        const DemoCity *ca = NULL, *cb = NULL;
        for (size_t i = 0; i < ncities; i++) {
            if (cities[i].id == corridors[c].a) ca = &cities[i];
            if (cities[i].id == corridors[c].b) cb = &cities[i];
        }
        if (!ca || !cb) continue;

        int pts = corridors[c].way_pts;
        int64_t prev = ca->id;
        for (int p = 1; p <= pts; p++) {
            int64_t cur;
            if (p == pts) {
                cur = cb->id;
            } else {
                double t = (double)p / (double)pts;
                /* Add gentle organic curve bend to exercise degree-2 tortuosity */
                double curve = sin(t * M_PI) * 0.04 * ((c % 2 == 0) ? 1.0 : -1.0);
                double lat = ca->lat + t * (cb->lat - ca->lat) + curve;
                double lon = ca->lon + t * (cb->lon - ca->lon) - curve * 0.6;
                cur = next_id++;
                sqlite3_reset(nstmt);
                sqlite3_bind_int64(nstmt, 1, cur);
                sqlite3_bind_double(nstmt, 2, lat);
                sqlite3_bind_double(nstmt, 3, lon);
                sqlite3_step(nstmt);
            }

            /* Two-way carriageway edges */
            for (int dir = 0; dir < 2; dir++) {
                sqlite3_reset(estmt);
                sqlite3_bind_int64(estmt, 1, dir == 0 ? prev : cur);
                sqlite3_bind_int64(estmt, 2, dir == 0 ? cur : prev);
                sqlite3_bind_int(estmt, 3, (int)corridors[c].rc);
                sqlite3_bind_double(estmt, 4, corridors[c].spd);
                sqlite3_bind_double(estmt, 5, corridors[c].smooth);
                sqlite3_bind_int(estmt, 6, corridors[c].flags);
                sqlite3_step(estmt);
            }
            prev = cur;
        }
    }

    /* Add 10 cul-de-sac residential stubs to verify 1-core peeling */
    for (int s = 0; s < 10; s++) {
        int64_t spur_id = next_id++;
        sqlite3_reset(nstmt);
        sqlite3_bind_int64(nstmt, 1, spur_id);
        sqlite3_bind_double(nstmt, 2, cities[s].lat + 0.015);
        sqlite3_bind_double(nstmt, 3, cities[s].lon - 0.015);
        sqlite3_step(nstmt);

        for (int dir = 0; dir < 2; dir++) {
            sqlite3_reset(estmt);
            sqlite3_bind_int64(estmt, 1, dir == 0 ? cities[s].id : spur_id);
            sqlite3_bind_int64(estmt, 2, dir == 0 ? spur_id : cities[s].id);
            sqlite3_bind_int(estmt, 3, (int)RC_RESIDENTIAL);
            sqlite3_bind_double(estmt, 4, 32.2);
            sqlite3_bind_double(estmt, 5, 0.70);
            sqlite3_bind_int(estmt, 6, 0);
            sqlite3_step(estmt);
        }
    }

    sqlite3_finalize(nstmt);
    sqlite3_finalize(estmt);
    sql_exec(db, "COMMIT;");
    LOG_INFO("Seeded UK-wide demo road network (~%lld raw nodes)", (long long)next_id - 1000 + (long long)ncities);
    return 0;
}

static int db_load_raw_graph(sqlite3 *db, VmrRawGraph *rg) {
    raw_graph_free(rg);
    raw_graph_init(rg);

    /* Read-path pragmas: large page cache + mmap the file + in-memory temp.
     * These cut the row-by-row load time substantially on multi-GB DBs. */
    sqlite3_exec(db, "PRAGMA cache_size = -262144;", NULL, NULL, NULL);   /* 256 MB */
    sqlite3_exec(db, "PRAGMA mmap_size = 1073741824;", NULL, NULL, NULL);  /* 1 GB */
    sqlite3_exec(db, "PRAGMA temp_store = MEMORY;", NULL, NULL, NULL);
    sqlite3_exec(db, "PRAGMA query_only = ON;", NULL, NULL, NULL);

    sqlite3_stmt *nstmt = NULL;
    if (sqlite3_prepare_v2(db, "SELECT id, lat, lon FROM osm_nodes;", -1, &nstmt, NULL) != SQLITE_OK) {
        LOG_ERR("Failed to query osm_nodes: %s", sqlite3_errmsg(db));
        return -1;
    }
    while (sqlite3_step(nstmt) == SQLITE_ROW) {
        int64_t id = sqlite3_column_int64(nstmt, 0);
        double lat = sqlite3_column_double(nstmt, 1);
        double lon = sqlite3_column_double(nstmt, 2);
        if (raw_graph_add_node(rg, id, lat, lon) != 0) {
            sqlite3_finalize(nstmt);
            return -1;
        }
    }
    sqlite3_finalize(nstmt);

    sqlite3_stmt *estmt = NULL;
    if (sqlite3_prepare_v2(
            db, "SELECT from_node, to_node, road_class, speed_kmh, smooth, flags FROM osm_edges;",
            -1, &estmt, NULL) != SQLITE_OK) {
        LOG_ERR("Failed to query osm_edges: %s", sqlite3_errmsg(db));
        return -1;
    }

    while (sqlite3_step(estmt) == SQLITE_ROW) {
        int64_t fid = sqlite3_column_int64(estmt, 0);
        int64_t tid = sqlite3_column_int64(estmt, 1);
        uint32_t u = idmap_get(&rg->osm_map, fid);
        uint32_t v = idmap_get(&rg->osm_map, tid);
        if (u == VMR_INVALID || v == VMR_INVALID || u == v) {
            continue;
        }
        int rc_int = sqlite3_column_int(estmt, 2);
        RoadClass rc = (rc_int >= 0 && rc_int <= RC_UNCLASSIFIED) ? (RoadClass)rc_int : RC_UNCLASSIFIED;
        float spd_kmh = (float)sqlite3_column_double(estmt, 3);
        if (spd_kmh < 5.0f) spd_kmh = 30.0f;
        float smooth = (float)sqlite3_column_double(estmt, 4);
        if (smooth <= 0.0f || smooth > 1.0f) smooth = 0.85f;
        uint8_t flags = (uint8_t)sqlite3_column_int(estmt, 5);

        float def_spd = 48.3f, f_class = 0.30f;
        road_class_defaults(rc, &def_spd, &f_class);
        float f_speed = spd_kmh / 112.65f;
        if (f_speed > 1.0f) f_speed = 1.0f;

        double dist_m = haversine_m(*(rg->node_lat + u), *(rg->node_lon + u),
                                    *(rg->node_lat + v), *(rg->node_lon + v));
        if (dist_m < 0.5) dist_m = 0.5;
        float speed_mps = spd_kmh / 3.6f;
        float time_s = (float)dist_m / speed_mps;

        if (raw_graph_add_edge(rg, u, v, (float)dist_m, time_s,
                               f_speed, f_class, smooth, (uint8_t)rc, flags) != 0) {
            sqlite3_finalize(estmt);
            return -1;
        }
    }
    sqlite3_finalize(estmt);

    sqlite3_stmt *tstmt = NULL;
    if (sqlite3_prepare_v2(db, "SELECT via_node,from_node,to_node,mode FROM turn_restrictions;",
                           -1, &tstmt, NULL) == SQLITE_OK) {
        while (sqlite3_step(tstmt) == SQLITE_ROW) {
            int64_t via = sqlite3_column_int64(tstmt, 0);
            int64_t from = sqlite3_column_int64(tstmt, 1);
            int64_t to = sqlite3_column_int64(tstmt, 2);
            int mode = sqlite3_column_int(tstmt, 3);
            uint32_t vi = idmap_get(&rg->osm_map, via);
            uint32_t fi = idmap_get(&rg->osm_map, from);
            uint32_t ti = idmap_get(&rg->osm_map, to);
            if (vi == VMR_INVALID || fi == VMR_INVALID || ti == VMR_INVALID || vi == fi || vi == ti) {
                continue;
            }
            if (rg->turn_count % 1024U == 0) {
                uint32_t ncap = rg->turn_count ? rg->turn_count * 2U : 1024U;
                uint32_t *tv = (uint32_t *)xrealloc(rg->turn_via, (size_t)ncap * sizeof(uint32_t));
                uint32_t *tf = (uint32_t *)xrealloc(rg->turn_from, (size_t)ncap * sizeof(uint32_t));
                uint32_t *tt = (uint32_t *)xrealloc(rg->turn_to, (size_t)ncap * sizeof(uint32_t));
                uint8_t  *tm = (uint8_t *)xrealloc(rg->turn_mode, (size_t)ncap);
                if (tv) rg->turn_via = tv;
                if (tf) rg->turn_from = tf;
                if (tt) rg->turn_to = tt;
                if (tm) rg->turn_mode = tm;
                if (!tv || !tf || !tt || !tm) {
                    sqlite3_finalize(tstmt);
                    return -1;
                }
            }
            *(rg->turn_via + rg->turn_count) = vi;
            *(rg->turn_from + rg->turn_count) = fi;
            *(rg->turn_to + rg->turn_count) = ti;
            *(rg->turn_mode + rg->turn_count) = (uint8_t)(mode ? 1 : 0);
            rg->turn_count++;
        }
        sqlite3_finalize(tstmt);
    }

    LOG_INFO("Loaded raw SQLite graph: %u nodes, %u directed edges, %u turn restrictions",
             rg->node_count, rg->edge_count, rg->turn_count);
    return 0;
}

/* ========================================================================= */
/* SECTION 8: SPATIAL SNAP & 8D VECTOR-CORRIDOR BIDIRECTIONAL PATHFINDER     */
/* ========================================================================= */

static uint32_t vmr_snap_coord(const VmrGraph *g, double lat, double lon, int prefer_non_dead) {
    if (g->node_count == 0) {
        return VMR_INVALID;
    }
    float qx = 0.0f, qy = 0.0f;
    uk_project_tm(lat, lon, &qx, &qy);

    float span_x = g->grid_max_x - g->grid_min_x;
    float span_y = g->grid_max_y - g->grid_min_y;
    if (span_x < 1.0f) span_x = 1.0f;
    if (span_y < 1.0f) span_y = 1.0f;

    int cx = (int)(((qx - g->grid_min_x) / span_x) * (float)(SPATIAL_GRID_DIM - 1));
    int cy = (int)(((qy - g->grid_min_y) / span_y) * (float)(SPATIAL_GRID_DIM - 1));
    if (cx < 0) cx = 0; else if (cx >= SPATIAL_GRID_DIM) cx = SPATIAL_GRID_DIM - 1;
    if (cy < 0) cy = 0; else if (cy >= SPATIAL_GRID_DIM) cy = SPATIAL_GRID_DIM - 1;

    uint32_t best = VMR_INVALID;
    float best_d2 = FLT_MAX;

    /* Two passes: prefer live intersections, then accept a peeled cul-de-sac
     * node if the whole neighbourhood is dead (otherwise an all-dead island or
     * stub component could never be used as an endpoint). */
for (int pass = 0; pass < 2 && best == VMR_INVALID; pass++) {
        int skip_dead = (pass == 0) ? prefer_non_dead : 0;

        /* Expand search rings around target grid cell. When a ring first yields
         * a candidate, scan one extra ring before stopping: a query near a cell
         * boundary must not snap to a node kilometres away in its own cell while
         * a closer node sits in the neighbouring cell. */
        int found_rad = -1;
        for (int rad = 0; rad <= 8; rad++) {
            if (found_rad >= 0 && rad > found_rad + 1) {
                break;
            }
            int y0 = cy - rad; if (y0 < 0) y0 = 0;
            int y1 = cy + rad; if (y1 >= SPATIAL_GRID_DIM) y1 = SPATIAL_GRID_DIM - 1;
            int x0 = cx - rad; if (x0 < 0) x0 = 0;
            int x1 = cx + rad; if (x1 >= SPATIAL_GRID_DIM) x1 = SPATIAL_GRID_DIM - 1;

            for (int gy = y0; gy <= y1; gy++) {
                for (int gx = x0; gx <= x1; gx++) {
                    uint32_t cell = (uint32_t)(gy * SPATIAL_GRID_DIM + gx);
                    uint32_t beg = *(g->grid_offsets + cell);
                    uint32_t end = *(g->grid_offsets + cell + 1U);
                    for (uint32_t i = beg; i < end; i++) {
                        uint32_t ni = *(g->grid_nodes + i);
                        if (skip_dead && (*(g->node_flags + ni) & NODE_DEAD_TREE)) {
                            continue;
                        }
                        float dx = *(g->node_x + ni) - qx;
                        float dy = *(g->node_y + ni) - qy;
                        float d2 = dx * dx + dy * dy;
                        if (d2 < best_d2) {
                            best_d2 = d2;
                            best = ni;
                        }
                    }
                }
            }
            if (best != VMR_INVALID && found_rad < 0) {
                found_rad = rad;
            }
        }
    }

    /* Fallback global scan if query coordinate lies far outside grid cell radius */
    if (best == VMR_INVALID) {
        for (uint32_t ni = 0; ni < g->node_count; ni++) {
            float dx = *(g->node_x + ni) - qx;
            float dy = *(g->node_y + ni) - qy;
            float d2 = dx * dx + dy * dy;
            if (d2 < best_d2) {
                best_d2 = d2;
                best = ni;
            }
        }
    }
    return best;
}

static VmrPreference vmr_preset_profile(const char *name) {
    VmrPreference p;
    int id = 0;
    if (name == NULL || c_stricmp(name, "balanced") == 0) {
        p = (VmrPreference){ 1.40f, 1.20f, 1.00f, 0.60f, 0.50f, 0, 0 };
        id = 0;
    } else if (c_stricmp(name, "fastest") == 0) {
        p = (VmrPreference){ 0.80f, 2.20f, 1.50f, 0.30f, 0.40f, 0, 0 };
        id = 1;
    } else if (c_stricmp(name, "direct") == 0) {
        p = (VmrPreference){ 2.60f, 0.50f, 0.30f, 1.20f, 0.30f, 0, 0 };
        id = 2;
    } else if (c_stricmp(name, "highway") == 0) {
        p = (VmrPreference){ 1.10f, 1.80f, 2.20f, 0.50f, 0.60f, 0, 0 };
        id = 3;
    } else if (c_stricmp(name, "local") == 0 || c_stricmp(name, "no-toll") == 0) {
        p = (VmrPreference){ 1.30f, 0.70f, 0.20f, 0.80f, 1.40f, 1, 0 };
        id = 4;
    } else {
        p = (VmrPreference){ 1.40f, 1.20f, 1.00f, 0.60f, 0.50f, 0, 0 };
        id = 0;
    }
    p.pref_id = id;
    return p;
}

/*
 * Computes 2-hop lookahead cosine similarity using contiguous SoA arrays.
 * Adds candidate edge vector (dx1, dy1) with each non-reversing continuation
 * edge (dx2, dy2) at `v` and measures dot product against target unit vector
 * (rx, ry). Runs strictly over L1-resident floats!
 */
static inline float two_hop_fwd_cos(const VmrGraph *g, uint32_t u, uint32_t v,
                                      uint32_t goal, float dx1, float dy1,
                                      float cos1, float rx, float ry) {
    if (v == goal) {
        return 1.0f;
    }
    uint32_t beg = *(g->fwd_offsets + v);
    uint32_t end = *(g->fwd_offsets + v + 1U);
    if (end - beg > 8U) {
        return cos1; /* Skip expensive fanout on giant hubs */
    }
    float best = cos1;
    int found = 0;
    const uint32_t *to_p = g->edge_to;
    const float *dx_p = g->edge_dx;
    const float *dy_p = g->edge_dy;

    for (uint32_t e2 = beg; e2 < end; e2++) {
        uint32_t w = *(to_p + e2);
        if (w == u) continue;
        if (w == goal) return 1.0f;
        float sx = dx1 + *(dx_p + e2);
        float sy = dy1 + *(dy_p + e2);
        float inv_n = 1.0f / (sqrtf(sx * sx + sy * sy) + EPS_F);
        float c2 = (sx * rx + sy * ry) * inv_n;
        if (!found || c2 > best) {
            best = c2;
            found = 1;
        }
    }
    return best;
}

static inline float two_hop_rev_cos(const VmrGraph *g, uint32_t w, uint32_t u,
                                      uint32_t start, float dx1, float dy1,
                                      float cos1, float rx, float ry) {
    if (u == start) {
        return 1.0f;
    }
    uint32_t beg = *(g->rev_offsets + u);
    uint32_t end = *(g->rev_offsets + u + 1U);
    if (end - beg > 8U) {
        return cos1;
    }
    float best = cos1;
    int found = 0;
    const uint32_t *src_p = g->rev_src_node;
    const uint32_t *eid_p = g->rev_fwd_eid;
    const float *dx_p = g->edge_dx;
    const float *dy_p = g->edge_dy;

    for (uint32_t r2 = beg; r2 < end; r2++) {
        uint32_t prev = *(src_p + r2);
        if (prev == w) continue;
        if (prev == start) return 1.0f;
        uint32_t eid2 = *(eid_p + r2);
        /* Backward traversal flips edge vector sign */
        float sx = (-dx1) - *(dx_p + eid2);
        float sy = (-dy1) - *(dy_p + eid2);
        float inv_n = 1.0f / (sqrtf(sx * sx + sy * sy) + EPS_F);
        float c2 = (sx * rx + sy * ry) * inv_n;
        if (!found || c2 > best) {
            best = c2;
            found = 1;
        }
    }
    return best;
}

/*
 * Cache-friendly two-hop lookup: the continuation sums and inverse norms are
 * precomputed per edge, so this is a contiguous dot product with no adjacency
 * scan and no runtime sqrt. `e` is the forward edge id, `ri` the reverse index.
 */
static inline float two_hop_fwd_cos_fast(const VmrGraph *g, uint32_t e,
                                          float cos1, float rx, float ry) {
    uint32_t c0 = *(g->th_off + e), c1 = *(g->th_off + e + 1U);
    if (c0 == c1) {
        return cos1;
    }
    float best = cos1;
    int found = 0;
    const float *sx = g->th_sx + c0, *sy = g->th_sy + c0, *iv = g->th_inv + c0;
    for (uint32_t c = c0; c < c1; c++) {
        uint32_t k = c - c0;
        float cc = (*(sx + k) * rx + *(sy + k) * ry) * *(iv + k);
        if (!found || cc > best) {
            best = cc;
            found = 1;
        }
    }
    return best;
}

static inline float two_hop_rev_cos_fast(const VmrGraph *g, uint32_t ri,
                                          float cos1, float rx, float ry) {
    uint32_t c0 = *(g->rth_off + ri), c1 = *(g->rth_off + ri + 1U);
    if (c0 == c1) {
        return cos1;
    }
    float best = cos1;
    int found = 0;
    const float *sx = g->rth_sx + c0, *sy = g->rth_sy + c0, *iv = g->rth_inv + c0;
    for (uint32_t c = c0; c < c1; c++) {
        uint32_t k = c - c0;
        float cc = (*(sx + k) * rx + *(sy + k) * ry) * *(iv + k);
        if (!found || cc > best) {
            best = cc;
            found = 1;
        }
    }
    return best;
}

/*
 * Computes the 8-Dimensional Vector Edge Cost:
 *   Spatial dims:  (1) unit dot product cos_1, (2) 2-hop cos_2, (3) momentum dot
 *   Feature dims:  (4) speed, (5) class, (6) straightness, (7) smoothness,
 *                  (8) cross-track corridor rejection relative to A -> B chord
 */
static inline float compute_attr_scale(const VmrGraph *g, uint32_t eid,
                                       const VmrPreference *p) {
    float dot = (p->w_speed    * *(g->edge_f_speed + eid)) +
                (p->w_class    * *(g->edge_f_class + eid)) +
                (p->w_straight * *(g->edge_f_straight + eid)) +
                (p->w_smooth   * *(g->edge_f_smooth + eid));
    return 1.0f / (1.0f + 0.22f * dot);
}

static inline float compute_8d_step_cost(const VmrGraph *g, uint32_t eid,
                                         const VmrPreference *p, SearchMode mode,
                                         float attr_scale,
                                         float cos1, float cos2, float mom_dot,
                                         float cross_ratio) {
    float base_time = *(g->edge_time_s + eid);
    if (p->avoid_tolls && (*(g->edge_flags + eid) & EDGE_TOLL)) {
        base_time *= 3.8f;
    }
    if (mode == MODE_DIJKSTRA || mode == MODE_ASTAR) {
        return base_time;
    }

    /* Spatial vector alignment (1-hop + 2-hop sum + directional momentum) */
    float spatial_align = (0.50f * cos1) + (0.30f * cos2) + (0.20f * mom_dot);
    if (spatial_align > 1.0f) spatial_align = 1.0f;
    if (spatial_align < -1.0f) spatial_align = -1.0f;

    float align_pen = 1.0f + (0.38f * p->w_align * (1.0f - spatial_align));

    /* Cross-track vector rejection penalty (keeps route in A<->B corridor) */
    if (cross_ratio > 1.5f) cross_ratio = 1.5f;
    float corridor_pen = 1.0f + (0.14f * p->w_align * cross_ratio);

    float cost = base_time * attr_scale * align_pen * corridor_pen;
    return (cost > 0.01f) ? cost : 0.01f;
}

/* Cheaper admissible lower bound: full cost with cos2=1 and corridor term = 1. */
static inline float compute_8d_cost_min(const VmrGraph *g, uint32_t eid,
                                         const VmrPreference *p, SearchMode mode,
                                         float attr_scale, float cos1, float mom_dot) {
    float base_time = *(g->edge_time_s + eid);
    if (p->avoid_tolls && (*(g->edge_flags + eid) & EDGE_TOLL)) {
        base_time *= 3.8f;
    }
    if (mode == MODE_DIJKSTRA || mode == MODE_ASTAR) {
        return base_time;
    }
    float spatial_align = 0.50f * cos1 + 0.30f + 0.20f * mom_dot;
    if (spatial_align > 1.0f) spatial_align = 1.0f;
    if (spatial_align < -1.0f) spatial_align = -1.0f;
    float align_pen = 1.0f + (0.38f * p->w_align * (1.0f - spatial_align));
    float cost = base_time * attr_scale * align_pen;
    return (cost > 0.01f) ? cost : 0.01f;
}

static inline float alt_heuristic_scalar(const AltNode32 *av, const float *gf,
                                         const float *gb, uint32_t L) {
    float best = 0.0f;
    for (uint32_t i = 0; i < L && i < 32; i++) {
        float d1 = *(gf + i) - av->fwd[i];
        float d2 = av->bwd[i] - *(gb + i);
        if (d1 > best) best = d1;
        if (d2 > best) best = d2;
    }
    return best;
}

#if defined(__AVX2__)
#include <immintrin.h>
static inline float alt_heuristic_avx2_8(const AltNode32 *av, __m256 gf0, __m256 gb0) {
    __m256 f0 = _mm256_loadu_ps(av->fwd);
    __m256 b0 = _mm256_loadu_ps(av->bwd);
    __m256 d1 = _mm256_sub_ps(gf0, f0);
    __m256 d2 = _mm256_sub_ps(b0, gb0);
    __m256 mx = _mm256_max_ps(_mm256_max_ps(d1, d2), _mm256_setzero_ps());
    __m128 hi = _mm256_extractf128_ps(mx, 1);
    __m128 lo = _mm256_castps256_ps128(mx);
    __m128 m = _mm_max_ps(lo, hi);
    m = _mm_max_ps(m, _mm_movehl_ps(m, m));
    m = _mm_max_ss(m, _mm_shuffle_ps(m, m, 1));
    return _mm_cvtss_f32(m);
}
static inline float alt_heuristic_avx2(const AltNode32 *av, __m256 gf0, __m256 gf1,
                                       __m256 gf2, __m256 gf3,
                                       __m256 gb0, __m256 gb1, __m256 gb2, __m256 gb3) {
    __m256 f0 = _mm256_loadu_ps(av->fwd);
    __m256 f1 = _mm256_loadu_ps(av->fwd + 8);
    __m256 f2 = _mm256_loadu_ps(av->fwd + 16);
    __m256 f3 = _mm256_loadu_ps(av->fwd + 24);
    __m256 b0 = _mm256_loadu_ps(av->bwd);
    __m256 b1 = _mm256_loadu_ps(av->bwd + 8);
    __m256 b2 = _mm256_loadu_ps(av->bwd + 16);
    __m256 b3 = _mm256_loadu_ps(av->bwd + 24);
    __m256 d1 = _mm256_max_ps(_mm256_max_ps(_mm256_sub_ps(gf0, f0), _mm256_sub_ps(gf1, f1)),
                              _mm256_max_ps(_mm256_sub_ps(gf2, f2), _mm256_sub_ps(gf3, f3)));
    __m256 d2 = _mm256_max_ps(_mm256_max_ps(_mm256_sub_ps(b0, gb0), _mm256_sub_ps(b1, gb1)),
                              _mm256_max_ps(_mm256_sub_ps(b2, gb2), _mm256_sub_ps(b3, gb3)));
    __m256 mx = _mm256_max_ps(_mm256_max_ps(d1, d2), _mm256_setzero_ps());
    __m128 hi = _mm256_extractf128_ps(mx, 1);
    __m128 lo = _mm256_castps256_ps128(mx);
    __m128 m = _mm_max_ps(lo, hi);
    m = _mm_max_ps(m, _mm_movehl_ps(m, m));
    m = _mm_max_ss(m, _mm_shuffle_ps(m, m, 1));
    return _mm_cvtss_f32(m);
}
#endif

#if defined(__AVX2__)
#define ALT_HEUR(g, u, gf, gb, gf0, gf1, gf2, gf3, gb0, gb1, gb2, gb3, L) alt_heuristic_avx2_8(&(g)->alt_nodes[u], (gf0), (gb0))
#else
#define ALT_HEUR(g, u, gf, gb, gf0, gf1, gf2, gf3, gb0, gb1, gb2, gb3, L) alt_heuristic_scalar(&(g)->alt_nodes[u], (gf), (gb), (L))
#endif

static inline int turn_forbidden(const VmrGraph *g, uint32_t via,
                                 uint32_t from, uint32_t to) {
    if (g->turn_count == 0 || g->turn_offsets == NULL) {
        return 0;
    }
    uint32_t beg = *(g->turn_offsets + via);
    uint32_t end = *(g->turn_offsets + via + 1U);
    for (uint32_t i = beg; i < end; i++) {
        if (*(g->turn_from + i) == from && *(g->turn_to + i) == to) {
            return (*(g->turn_mode + i) == 0);
        }
    }
    return 0;
}

static void vmr_compute_result_telemetry(const VmrGraph *g, VmrResult *res,
                                         const VmrPreference *pref) {
    if (!res->found || res->hop_count == 0) {
        return;
    }
    double tot_dist = 0.0;
    double tot_time = 0.0;
    double sum_cos_w = 0.0;
    double sum_8d_w  = 0.0;
    double net_dx = 0.0;
    double net_dy = 0.0;

    float gx = *(g->node_x + res->goal_idx);
    float gy = *(g->node_y + res->goal_idx);
    float sx = *(g->node_x + res->start_idx);
    float sy = *(g->node_y + res->start_idx);
    double chord_dx = (double)(gx - sx);
    double chord_dy = (double)(gy - sy);
    double chord_len = sqrt(chord_dx * chord_dx + chord_dy * chord_dy);
    res->direct_distance_m = chord_len;

    const uint32_t *ep = res->edges;
    const uint32_t *eend = res->edges + res->hop_count;

    while (ep < eend) {
        uint32_t eid = *ep;
        uint32_t u = *(g->edge_from + eid);
        float len = *(g->edge_length_m + eid);
        float tm  = *(g->edge_time_s + eid);
        float dx  = *(g->edge_dx + eid);
        float dy  = *(g->edge_dy + eid);
        float ux  = *(g->edge_ux + eid);
        float uy  = *(g->edge_uy + eid);

        tot_dist += (double)len;
        tot_time += (double)tm;
        net_dx   += (double)dx;
        net_dy   += (double)dy;

        float rx = gx - *(g->node_x + u);
        float ry = gy - *(g->node_y + u);
        float rnorm = sqrtf(rx * rx + ry * ry);
        float local_cos = (rnorm > EPS_F) ? ((ux * rx + uy * ry) / rnorm) : 1.0f;
        sum_cos_w += (double)local_cos * (double)len;

        float attr_norm = (pref->w_speed * *(g->edge_f_speed + eid) +
                           pref->w_class * *(g->edge_f_class + eid) +
                           pref->w_straight * *(g->edge_f_straight + eid) +
                           pref->w_smooth * *(g->edge_f_smooth + eid)) /
                          (pref->w_speed + pref->w_class + pref->w_straight + pref->w_smooth + EPS_F);
        float align8 = 0.55f * local_cos + 0.45f * attr_norm;
        sum_8d_w += (double)align8 * (double)len;
        ep++;
    }

    double net_len = sqrt(net_dx * net_dx + net_dy * net_dy);
    double cum_cos = (chord_len > 1e-6 && net_len > 1e-6)
                     ? ((net_dx * chord_dx + net_dy * chord_dy) / (net_len * chord_len))
                     : 1.0;

    res->total_distance_m = tot_dist;
    res->est_time_s       = tot_time;
    res->detour_ratio     = (chord_len > 1e-3) ? (tot_dist / chord_len) : 1.0;
    res->mean_cos_sim     = (tot_dist > 1e-3) ? (sum_cos_w / tot_dist) : 1.0;
    res->cum_vec_cos_sim  = cum_cos;
    res->mean_8d_alignment= (tot_dist > 1e-3) ? (sum_8d_w / tot_dist) : 1.0;
}

/* Reusable search workspace: avoids ~14 allocations and a ~65 MB memset per
 * query.  Generation stamps make visited/closed tests O(1) without clearing.
 * The CLI is single-threaded; this is intentionally a process-global pool. */
typedef struct {
    float g;
    uint32_t p, pe, st;
    float mx, my;
} BidirFwdState;

typedef struct {
    float g;
    uint32_t p, pe, st;
    float mx, my;
} BidirBwdState;

_Static_assert(sizeof(BidirFwdState) == 24, "BidirFwdState must be 24 bytes");
_Static_assert(sizeof(BidirBwdState) == 24, "BidirBwdState must be 24 bytes");

typedef struct {
    size_t cap;
    int ready;
    uint32_t epoch;
    BidirFwdState *sf;
    BidirBwdState *sb;
    /* Legacy arrays retained only for the dormant CH helpers. */
    float *g_fwd, *g_bwd;
    float *mx_fwd, *my_fwd, *mx_bwd, *my_bwd;
    uint32_t *p_fwd, *p_bwd;
    uint32_t *pe_fwd, *pe_bwd;
    uint32_t *st_fwd, *cl_fwd, *st_bwd, *cl_bwd;
    MinHeap pq_fwd, pq_bwd;
} VmrSearchWS;

static VmrSearchWS g_ws;
static int g_ws_atexit = 0;

static void vmr_ws_free(void) {
    xfree(g_ws.g_fwd); xfree(g_ws.g_bwd);
    xfree(g_ws.mx_fwd); xfree(g_ws.my_fwd);
    xfree(g_ws.mx_bwd); xfree(g_ws.my_bwd);
    xfree(g_ws.p_fwd); xfree(g_ws.p_bwd);
    xfree(g_ws.pe_fwd); xfree(g_ws.pe_bwd);
    xfree(g_ws.st_fwd); xfree(g_ws.cl_fwd);
    xfree(g_ws.st_bwd); xfree(g_ws.cl_bwd);
    heap_free(&g_ws.pq_fwd);
    heap_free(&g_ws.pq_bwd);
    xfree_aligned(g_ws.sf);
    xfree_aligned(g_ws.sb);
    memset(&g_ws, 0, sizeof(g_ws));
}

static int vmr_ws_reserve(size_t n) {
    if (!g_ws_atexit) {
        atexit(vmr_ws_free);
        g_ws_atexit = 1;
    }
    if (g_ws.ready && g_ws.cap >= n) {
        return 0;
    }
    xfree(g_ws.g_fwd); xfree(g_ws.g_bwd);
    xfree(g_ws.mx_fwd); xfree(g_ws.my_fwd);
    xfree(g_ws.mx_bwd); xfree(g_ws.my_bwd);
    xfree(g_ws.p_fwd); xfree(g_ws.p_bwd);
    xfree(g_ws.pe_fwd); xfree(g_ws.pe_bwd);
    xfree(g_ws.st_fwd); xfree(g_ws.cl_fwd);
    xfree(g_ws.st_bwd); xfree(g_ws.cl_bwd);
    heap_free(&g_ws.pq_fwd);
    heap_free(&g_ws.pq_bwd);
    xfree_aligned(g_ws.sf);
    xfree_aligned(g_ws.sb);
    memset(&g_ws, 0, sizeof(g_ws));
    if (n == 0) {
        n = 1;
    }
    g_ws.cap = n;
    g_ws.sf = (BidirFwdState *)xaligned_alloc(64, n * sizeof(BidirFwdState));
    g_ws.sb = (BidirBwdState *)xaligned_alloc(64, n * sizeof(BidirBwdState));
    if (!g_ws.sf || !g_ws.sb || heap_init(&g_ws.pq_fwd, 256) != 0 || heap_init(&g_ws.pq_bwd, 256) != 0) {
        return -1;
    }
    g_ws.ready = 1;
    return 0;
}

/* Corridor around the scalar-CH path: a uniform grid of path nodes. */
static float g_corr_x1[128], g_corr_y1[128], g_corr_x2[128], g_corr_y2[128];
static uint32_t g_corr_n = 0;
static float g_corr_pad2 = 0.0f;
static int g_corr_disabled = 0;

static void vmr_corridor_free(void) {
    g_corr_n = 0;
    g_corr_pad2 = 0.0f;
}

static inline int in_corridor(const VmrGraph *g, uint32_t u) {
    if (g_corr_n == 0) {
        return 1;
    }
    float px = *(g->node_x + u);
    float py = *(g->node_y + u);
    for (uint32_t i = 0; i < g_corr_n; i++) {
        float x1 = *(g_corr_x1 + i), y1 = *(g_corr_y1 + i);
        float x2 = *(g_corr_x2 + i), y2 = *(g_corr_y2 + i);
        float dx = x2 - x1, dy = y2 - y1;
        float l2 = dx * dx + dy * dy;
        float t = (l2 > EPS_F) ? (((px - x1) * dx + (py - y1) * dy) / l2) : 0.0f;
        if (t < 0.0f) t = 0.0f;
        if (t > 1.0f) t = 1.0f;
        float qx = x1 + t * dx, qy = y1 + t * dy;
        float ex = px - qx, ey = py - qy;
        if (ex * ex + ey * ey <= g_corr_pad2) {
            return 1;
        }
    }
    return 0;
}

static __attribute__((unused)) int ch_nodes_push(uint32_t **nodes, uint32_t *n, uint32_t *cap, uint32_t v) {
    if (*n == *cap) {
        uint32_t ncap = *cap * 2U;
        uint32_t *nn = (uint32_t *)xrealloc(*nodes, (size_t)ncap * sizeof(uint32_t));
        if (!nn) {
            return -1;
        }
        *nodes = nn;
        *cap = ncap;
    }
    (*nodes)[(*n)++] = v;
    return 0;
}

static __attribute__((unused)) int ch_unpack_edge(const VmrGraph *g, uint32_t eid,
                          uint32_t **nodes, uint32_t *n, uint32_t *cap) {
    if (!g->ch_via || eid >= g->ch_edge_count) {
        LOG_ERR("CH unpack: bad edge id %u (edges %u)", eid, g->ch_edge_count);
        return -1;
    }
    uint32_t *stack = (uint32_t *)xmalloc(256 * sizeof(uint32_t));
    size_t sp = 0, scap = 256;
    if (!stack) {
        return -1;
    }
    *(stack + sp++) = eid;
    int rc = 0;
    while (sp > 0) {
        uint32_t e = *(stack + --sp);
        uint32_t via = *(g->ch_via + e);
        if (via == VMR_INVALID) {
            if (ch_nodes_push(nodes, n, cap, *(g->ch_to + e)) != 0) {
                rc = -1;
                break;
            }
            continue;
        }
        if (sp + 2 > scap) {
            size_t ncap = scap * 2U;
            uint32_t *ns = (uint32_t *)xrealloc(stack, ncap * sizeof(uint32_t));
            if (!ns) {
                rc = -1;
                break;
            }
            stack = ns;
            scap = ncap;
        }
        *(stack + sp++) = *(g->ch_child2 + e);
        *(stack + sp++) = *(g->ch_child1 + e);
    }
    xfree(stack);
    return rc;
}

static __attribute__((unused)) int ch_collect_edge(const VmrGraph *g, uint32_t eid,
                           uint32_t **nodes, uint32_t *n, uint32_t *cap,
                           uint32_t *vis) {
    if (eid >= g->ch_edge_count) {
        return 0;
    }
    uint32_t *stack = (uint32_t *)xmalloc(256 * sizeof(uint32_t));
    if (!stack) {
        return -1;
    }
    size_t sp = 0, scap = 256;
    *(stack + sp++) = eid;
    int rc = 0;
    while (sp > 0) {
        uint32_t e = *(stack + --sp);
        if (e >= g->ch_edge_count) {
            continue;
        }
        uint32_t w = e >> 5;
        uint32_t b = e & 31U;
        if ((*(vis + w) & (1U << b)) != 0) {
            continue;
        }
        *(vis + w) |= (1U << b);
        uint32_t via = *(g->ch_via + e);
        if (via == VMR_INVALID) {
            if (ch_nodes_push(nodes, n, cap, *(g->ch_to + e)) != 0) {
                rc = -1;
                break;
            }
            continue;
        }
        if (sp + 2 > scap) {
            size_t ncap = scap * 2U;
            uint32_t *ns = (uint32_t *)xrealloc(stack, ncap * sizeof(uint32_t));
            if (!ns) {
                rc = -1;
                break;
            }
            stack = ns;
            scap = ncap;
        }
        *(stack + sp++) = *(g->ch_child2 + e);
        *(stack + sp++) = *(g->ch_child1 + e);
    }
    xfree(stack);
    return rc;
}

/* Scalar travel-time CH search: returns the coarse path nodes (shortcut
 * endpoints) or NULL.  Used only to build the corridor. */
static __attribute__((unused)) uint32_t *vmr_ch_path(const VmrGraph *g, uint32_t start, uint32_t goal,
                             uint32_t *out_n) {
    *out_n = 0;
    if (!g->ch_built || start >= g->node_count || goal >= g->node_count) {
        return NULL;
    }
    uint32_t n = g->node_count;
    if (vmr_ws_reserve(n) != 0) {
        return NULL;
    }
    VmrSearchWS *ws = &g_ws;
    if (ws->epoch == UINT32_MAX) {
        memset(ws->st_fwd, 0, (size_t)n * sizeof(uint32_t));
        memset(ws->cl_fwd, 0, (size_t)n * sizeof(uint32_t));
        memset(ws->st_bwd, 0, (size_t)n * sizeof(uint32_t));
        memset(ws->cl_bwd, 0, (size_t)n * sizeof(uint32_t));
        ws->epoch = 0;
    }
    uint32_t epoch = ++ws->epoch;
    ws->pq_fwd.size = 0;
    ws->pq_bwd.size = 0;
    float *gf = ws->g_fwd, *gb = ws->g_bwd;
    uint32_t *pf = ws->p_fwd, *pb = ws->p_bwd;
    uint32_t *sf = ws->st_fwd, *sb = ws->st_bwd;
    uint32_t *cf = ws->cl_fwd, *cb = ws->cl_bwd;
    MinHeap *qf = &ws->pq_fwd, *qb = &ws->pq_bwd;

    *(sf + start) = epoch; *(gf + start) = 0.0f; *(pf + start) = VMR_INVALID;
    *(sb + goal) = epoch;  *(gb + goal) = 0.0f;  *(pb + goal) = VMR_INVALID;
    heap_push(qf, start, 0.0f);
    heap_push(qb, goal, 0.0f);

    uint32_t meet = VMR_INVALID;
    float best = FLT_MAX;
    int turn = 0;
    while (qf->size > 0 && qb->size > 0) {
        if (meet != VMR_INVALID) {
            float tf = qf->buf->prio + qb->buf->prio;
            if (tf >= best) {
                break;
            }
        }
        turn ^= 1;
        if (turn) {
            HeapEntry top;
            if (!heap_pop(qf, &top)) break;
            uint32_t u = top.node;
            if (*(cf + u) == epoch) continue;
            *(cf + u) = epoch;
            if (*(sb + u) == epoch && *(gb + u) < FLT_MAX) {
                float tot = *(gf + u) + *(gb + u);
                if (tot < best) { best = tot; meet = u; }
            }
            uint32_t b = *(g->ch_fwd_off + u), e = *(g->ch_fwd_off + u + 1U);
            for (; b < e; b++) {
                uint32_t v = *(g->ch_to + b);
                if (*(cf + v) == epoch) continue;
                if (!(*(g->ch_level + u) < *(g->ch_level + v))) continue;
                float ng = *(gf + u) + *(g->ch_time + b);
                if (*(sf + v) != epoch || ng < *(gf + v)) {
                    *(sf + v) = epoch;
                    *(gf + v) = ng;
                    *(pf + v) = u;
                    *(ws->pe_fwd + v) = b;
                    heap_push(qf, v, ng);
                }
            }
        } else {
            HeapEntry top;
            if (!heap_pop(qb, &top)) break;
            uint32_t w = top.node;
            if (*(cb + w) == epoch) continue;
            *(cb + w) = epoch;
            if (*(sf + w) == epoch && *(gf + w) < FLT_MAX) {
                float tot = *(gf + w) + *(gb + w);
                if (tot < best) { best = tot; meet = w; }
            }
            uint32_t b = *(g->ch_rev_off + w), e = *(g->ch_rev_off + w + 1U);
            for (; b < e; b++) {
                uint32_t u = *(g->ch_rev_src + b);
                uint32_t eid = *(g->ch_rev_eid + b);
                if (*(cb + u) == epoch) continue;
                if (!(*(g->ch_level + w) < *(g->ch_level + u))) continue;
                float ng = *(gb + w) + *(g->ch_time + eid);
                if (*(sb + u) != epoch || ng < *(gb + u)) {
                    *(sb + u) = epoch;
                    *(gb + u) = ng;
                    *(pb + u) = w;
                    *(ws->pe_bwd + u) = eid;
                    heap_push(qb, u, ng);
                }
            }
        }
    }
    if (meet == VMR_INVALID) {
        LOG_INFO("CH path: no meeting node");
        return NULL;
    }
    LOG_INFO("CH path: meet=%u best=%.2f", meet, best);
    uint32_t cur = meet;
    for (uint32_t k = 0; k < n && cur != start; k++) {
        cur = *(pf + cur);
        if (cur == VMR_INVALID) return NULL;
    }
    if (cur != start) {
        LOG_INFO("CH path: forward chain broken");
        return NULL;
    }
    LOG_INFO("CH path: forward chain ok");
    cur = meet;
    for (uint32_t k = 0; k < n && cur != goal; k++) {
        cur = *(pb + cur);
        if (cur == VMR_INVALID) return NULL;
    }
    if (cur != goal) {
        LOG_INFO("CH path: backward chain broken");
        return NULL;
    }
    LOG_INFO("CH path: backward chain ok");

    uint32_t fh = 0;
    cur = meet;
    while (cur != start && fh < n) {
        cur = *(pf + cur);
        if (cur == VMR_INVALID) break;
        fh++;
    }
    if (cur != start) {
        return NULL;
    }
    uint32_t bh = 0;
    cur = meet;
    while (cur != goal && bh < n) {
        cur = *(pb + cur);
        if (cur == VMR_INVALID) break;
        bh++;
    }
    if (cur != goal) {
        return NULL;
    }
    uint32_t total = fh + bh;
    uint32_t *nodes = (uint32_t *)xmalloc((size_t)(total + 1U) * sizeof(uint32_t));
    if (!nodes) {
        return NULL;
    }
    cur = meet;
    for (uint32_t k = fh;; k--) {
        *(nodes + k) = cur;
        if (k == 0) break;
        cur = *(pf + cur);
    }
    if (bh > 0) {
        cur = meet;
        for (uint32_t k = 0; k < bh; k++) {
            cur = *(pb + cur);
            *(nodes + fh + 1U + k) = cur;
        }
    }
    *out_n = total + 1U;
    return nodes;
}

/* Exact travel-time distances in the scalar CH graph.  These are lower bounds
 * for the 8D vector cost (cost >= time * attr_scale_min) and make the vector A*
 * goal-directed without changing its cost model. */
static __attribute__((unused)) void ch_distance_all(const VmrGraph *g, uint32_t source, float *dist, int reverse) {
    uint32_t n = g->node_count;
    if (!g->ch_built || !dist || source >= n) {
        return;
    }
    uint8_t *closed = (uint8_t *)xcalloc(n, 1);
    MinHeap pq;
    if (!closed || heap_init(&pq, 1024) != 0) {
        xfree(closed);
        return;
    }
    for (uint32_t i = 0; i < n; i++) {
        *(dist + i) = FLT_MAX;
    }
    *(dist + source) = 0.0f;
    heap_push(&pq, source, 0.0f);
    HeapEntry top;
    while (heap_pop(&pq, &top)) {
        uint32_t u = top.node;
        if (*(closed + u)) {
            continue;
        }
        *(closed + u) = 1;
        float du = *(dist + u);
        if (!reverse) {
            uint32_t b = *(g->ch_fwd_off + u), e = *(g->ch_fwd_off + u + 1U);
            for (; b < e; b++) {
                uint32_t v = *(g->ch_to + b);
                if (*(closed + v)) continue;
                if (!(*(g->ch_level + u) < *(g->ch_level + v))) continue;
                float nd = du + *(g->ch_time + b);
                if (nd < *(dist + v)) {
                    *(dist + v) = nd;
                    heap_push(&pq, v, nd);
                }
            }
        } else {
            uint32_t b = *(g->ch_rev_off + u), e = *(g->ch_rev_off + u + 1U);
            for (; b < e; b++) {
                uint32_t v = *(g->ch_rev_src + b);
                uint32_t eid = *(g->ch_rev_eid + b);
                if (*(closed + v)) continue;
                if (!(*(g->ch_level + v) < *(g->ch_level + u))) continue;
                float nd = du + *(g->ch_time + eid);
                if (nd < *(dist + v)) {
                    *(dist + v) = nd;
                    heap_push(&pq, v, nd);
                }
            }
        }
    }
    xfree(closed);
    heap_free(&pq);
}

#if defined(__AVX2__)
/*
 * AVX2 block expansion: four outgoing edges per block, cost computed 4-wide.
 * The two-hop cosines are resolved per lane from the precomputed continuation
 * table (sum vector + inverse norm), so the inner loop is a dot product.
 */
static inline uint32_t vmr_expand_fwd_block(
    const VmrGraph *g, VmrSearchWS *ws, MinHeap *pq, uint32_t epoch,
    uint32_t u, uint32_t goal, float gu, float umx, float umy,
    float urx, float ury, float gx, float gy, float sx, float sy,
    float chord_dx, float chord_dy, float inv_chord,
    const VmrPreference *pref, uint32_t pref_id, int use_cross,
    int use_alt, const float *alt_gf, const float *alt_gb, float alt_scale,
    uint32_t L, float v_heur, float *best_total, uint32_t *meet_node,
    int track_meet, int skip_dead) {
    uint32_t ebase = *(g->fwd_offsets + u);
    uint32_t end = *(g->fwd_offsets + u + 1U);
    if (ebase >= end || !g->edge_attr_scale[pref_id]) {
        return ebase;
    }
    const float ux_pos = *(g->node_x + u), uy_pos = *(g->node_y + u);

    __m128 urxv = _mm_set1_ps(urx), uryv = _mm_set1_ps(ury);
    __m128 umxv = _mm_set1_ps(umx), umyv = _mm_set1_ps(umy);
    __m128 sxf = _mm_set1_ps(sx), syf = _mm_set1_ps(sy);
    __m128 cdx = _mm_set1_ps(chord_dx), cdy = _mm_set1_ps(chord_dy);
    __m128 invc = _mm_set1_ps(inv_chord);
    __m128 wa = _mm_set1_ps(0.38f * pref->w_align);
    __m128 wc = _mm_set1_ps(0.14f * pref->w_align);
    __m128 zero = _mm_setzero_ps();
    __m128 one = _mm_set1_ps(1.0f);
    __m128 negone = _mm_set1_ps(-1.0f);
    __m128 half = _mm_set1_ps(0.5f), t3 = _mm_set1_ps(0.3f), t2 = _mm_set1_ps(0.2f);
    __m128 tollv = _mm_set1_ps(3.8f);
    __m128 mincost = _mm_set1_ps(0.01f);
    __m128 maxcross = _mm_set1_ps(1.5f);
    const float *th_sx = g->th_sx, *th_sy = g->th_sy, *th_inv = g->th_inv;
    const uint32_t *th_off = g->th_off;
    uint32_t e = ebase;

    for (; e + 4U <= end; e += 4U) {
        uint32_t tov[4];
        for (int l = 0; l < 4; l++) {
            tov[l] = *(g->edge_to + e + (uint32_t)l);
        }
        __m128 ux = _mm_loadu_ps(g->edge_ux + e), uy = _mm_loadu_ps(g->edge_uy + e);
        __m128 vxs = _mm_add_ps(_mm_set1_ps(ux_pos), _mm_loadu_ps(g->edge_dx + e));
        __m128 vys = _mm_add_ps(_mm_set1_ps(uy_pos), _mm_loadu_ps(g->edge_dy + e));
        __m128 tim = _mm_loadu_ps(g->edge_time_s + e);
        __m128 attr = _mm_loadu_ps(g->edge_attr_scale[pref_id] + e);
        __m128 cos1 = _mm_add_ps(_mm_mul_ps(ux, urxv), _mm_mul_ps(uy, uryv));
        __m128 mom = _mm_add_ps(_mm_mul_ps(umxv, ux), _mm_mul_ps(umyv, uy));

        __m128 cross = zero;
        if (use_cross) {
            __m128 cr = _mm_sub_ps(_mm_mul_ps(cdx, _mm_sub_ps(vys, syf)),
                                   _mm_mul_ps(cdy, _mm_sub_ps(vxs, sxf)));
            cr = _mm_andnot_ps(_mm_set1_ps(-0.0f), cr);
            cr = _mm_min_ps(_mm_mul_ps(_mm_mul_ps(cr, invc), invc), maxcross);
            cross = cr;
        }

        float cos1v[4];
        _mm_storeu_ps(cos1v, cos1);
        float c2v[4];
        for (int l = 0; l < 4; l++) {
            uint32_t v = tov[l];
            if (v == VMR_INVALID) {
                c2v[l] = cos1v[l];
                continue;
            }
            if (v == goal) {
                c2v[l] = 1.0f;
                continue;
            }
            uint32_t ei = e + (uint32_t)l;
            uint32_t c0 = *(th_off + ei), c1 = *(th_off + ei + 1U);
            float best = cos1v[l];
            int found = 0;
            for (uint32_t c = c0; c < c1; c++) {
                float cc = (*(th_sx + c) * urx + *(th_sy + c) * ury) * *(th_inv + c);
                if (!found || cc > best) {
                    best = cc;
                    found = 1;
                }
            }
            c2v[l] = found ? best : cos1v[l];
        }
        __m128 c2 = _mm_loadu_ps(c2v);

        __m128 spatial = _mm_add_ps(_mm_add_ps(_mm_mul_ps(half, cos1),
                                               _mm_mul_ps(t3, c2)),
                                    _mm_mul_ps(t2, mom));
        spatial = _mm_max_ps(_mm_min_ps(spatial, one), negone);
        __m128 align = _mm_add_ps(one, _mm_mul_ps(wa, _mm_sub_ps(one, spatial)));
        __m128 corr = _mm_add_ps(one, _mm_mul_ps(wc, cross));

        __m128 base;
        if (pref->avoid_tolls) {
            uint32_t f4;
            memcpy(&f4, g->edge_flags + e, sizeof(f4));
            __m128i fl = _mm_cvtepu8_epi32(_mm_cvtsi32_si128((int)f4));
            __m128i tm = _mm_and_si128(fl, _mm_set1_epi32(EDGE_TOLL));
            __m128 mask = _mm_castsi128_ps(_mm_cmpeq_epi32(tm, _mm_set1_epi32(EDGE_TOLL)));
            base = _mm_blendv_ps(tim, _mm_mul_ps(tim, tollv), mask);
        } else {
            base = tim;
        }

        __m128 cost = _mm_max_ps(_mm_mul_ps(_mm_mul_ps(base, attr),
                                            _mm_mul_ps(align, corr)), mincost);
        float costv[4];
        _mm_storeu_ps(costv, cost);

        __m128 spatial_min = _mm_add_ps(_mm_add_ps(_mm_mul_ps(half, cos1),
                                                    _mm_mul_ps(t3, one)),
                                         _mm_mul_ps(t2, mom));
        spatial_min = _mm_max_ps(_mm_min_ps(spatial_min, one), negone);
        __m128 align_min = _mm_add_ps(one, _mm_mul_ps(wa, _mm_sub_ps(one, spatial_min)));
        __m128 cost_min = _mm_max_ps(_mm_mul_ps(_mm_mul_ps(base, attr), align_min), mincost);
        float cminv[4];
        _mm_storeu_ps(cminv, cost_min);

        uint32_t p_u = ws->sf[u].p;
        int turn_res = ((*(g->node_flags + u) & NODE_HAS_TURN_RES) && p_u != VMR_INVALID);
        float gu_ = gu;
        for (int l = 0; l < 4; l++) {
            uint32_t v = tov[l];
            if (v == VMR_INVALID) {
                continue;
            }
            if (skip_dead && (*(g->edge_flags + e + (uint32_t)l) & EDGE_DEAD_TARGET)) {
                continue;
            }
            if (ws->sf[v].st == epoch + 1U) {
                continue;
            }
            if (turn_res && turn_forbidden(g, u, p_u, v)) {
                continue;
            }
            float h_v;
            if (use_alt) {
                h_v = alt_heuristic_scalar(&g->alt_nodes[v], alt_gf, alt_gb, L) * alt_scale;
            } else {
                float rvx = gx - *(g->node_x + v);
                float rvy = gy - *(g->node_y + v);
                h_v = vmr_sqrt_fast(rvx * rvx + rvy * rvy) * v_heur;
            }
            float ng_min = gu_ + cminv[l];
            if (ws->sf[v].st == epoch && ng_min >= ws->sf[v].g) {
                continue;
            }
            if (ng_min + h_v >= *best_total) {
                continue;
            }
            float ng = gu_ + costv[l];
            if (ng + h_v >= *best_total) {
                continue;
            }
            if (ws->sf[v].st != epoch || ng < ws->sf[v].g) {
                ws->sf[v].st = epoch;
                ws->sf[v].g = ng;
                ws->sf[v].p = u;
                ws->sf[v].pe = e + (uint32_t)l;
                float eux = *(g->edge_ux + e + (uint32_t)l);
                float euy = *(g->edge_uy + e + (uint32_t)l);
                float nmx = 0.65f * umx + 0.35f * eux;
                float nmy = 0.65f * umy + 0.35f * euy;
                float mn = vmr_inv_sqrt(nmx * nmx + nmy * nmy + EPS_F);
                ws->sf[v].mx = nmx * mn;
                ws->sf[v].my = nmy * mn;

                if (track_meet && ws->sb[v].st == epoch &&
                    ws->sb[v].g < FLT_MAX) {
                    float tot = ng + ws->sb[v].g;
                    if (tot < *best_total) {
                        *best_total = tot;
                        *meet_node = v;
                    }
                }
                if (ng + h_v < *best_total) {
                    heap_push(pq, v, ng + h_v);
                }
            }
        }
    }
    return e;
}
#endif

static inline uint32_t vmr_expand_bwd_block(
    const VmrGraph *g, VmrSearchWS *ws, MinHeap *pq, uint32_t epoch,
    uint32_t w, uint32_t start, float gb, float wmx, float wmy,
    float urx, float ury, float sx, float sy,
    float chord_dx, float chord_dy, float inv_chord,
    const VmrPreference *pref, uint32_t pref_id, int use_cross,
    int use_alt, const float *alt_sf, const float *alt_sb, float alt_scale,
    uint32_t L, float v_heur, float *best_total, uint32_t *meet_node,
    int track_meet, int skip_dead) {
    uint32_t bb0 = *(g->rev_block_off + w);
    uint32_t bb1 = *(g->rev_block_off + w + 1U);
    uint32_t rbase = *(g->rev_offsets + w);
    if (bb0 == bb1) {
        return rbase;
    }
    const EdgeBlock4 *blocks = g->rev_blocks;

    __m128 urxv = _mm_set1_ps(urx), uryv = _mm_set1_ps(ury);
    __m128 umxv = _mm_set1_ps(wmx), umyv = _mm_set1_ps(wmy);
    __m128 sxf = _mm_set1_ps(sx), syf = _mm_set1_ps(sy);
    __m128 cdx = _mm_set1_ps(chord_dx), cdy = _mm_set1_ps(chord_dy);
    __m128 invc = _mm_set1_ps(inv_chord);
    __m128 wa = _mm_set1_ps(0.38f * pref->w_align);
    __m128 wc = _mm_set1_ps(0.14f * pref->w_align);
    __m128 zero = _mm_setzero_ps();
    __m128 one = _mm_set1_ps(1.0f);
    __m128 negone = _mm_set1_ps(-1.0f);
    __m128 half = _mm_set1_ps(0.5f), t3 = _mm_set1_ps(0.3f), t2 = _mm_set1_ps(0.2f);
    __m128 tollv = _mm_set1_ps(3.8f);
    __m128 mincost = _mm_set1_ps(0.01f);
    __m128 maxcross = _mm_set1_ps(1.5f);
    const float *th_sx = g->rth_sx, *th_sy = g->rth_sy, *th_inv = g->rth_inv;
    const uint32_t *th_off = g->rth_off;

    for (uint32_t bb = bb0; bb < bb1; bb++) {
        const EdgeBlock4 *B = blocks + bb;
        __m128 ux = _mm_load_ps(B->ux), uy = _mm_load_ps(B->uy);
        __m128 vxs = _mm_load_ps(B->vx), vys = _mm_load_ps(B->vy);
        __m128 tim = _mm_load_ps(B->time);
        __m128 attr = _mm_load_ps(B->attr[pref_id]);
        __m128 cos1 = _mm_add_ps(_mm_mul_ps(ux, urxv), _mm_mul_ps(uy, uryv));
        __m128 mom = _mm_add_ps(_mm_mul_ps(umxv, ux), _mm_mul_ps(umyv, uy));

        __m128 cross = zero;
        if (use_cross) {
            __m128 cr = _mm_sub_ps(_mm_mul_ps(cdx, _mm_sub_ps(vys, syf)),
                                   _mm_mul_ps(cdy, _mm_sub_ps(vxs, sxf)));
            cr = _mm_andnot_ps(_mm_set1_ps(-0.0f), cr);
            cr = _mm_min_ps(_mm_mul_ps(_mm_mul_ps(cr, invc), invc), maxcross);
            cross = cr;
        }

        float cos1v[4];
        _mm_storeu_ps(cos1v, cos1);
        float c2v[4];
        uint32_t r0 = rbase + (bb - bb0) * 4U;
        for (int l = 0; l < 4; l++) {
            uint32_t v = B->to[l];
            if (v == VMR_INVALID) {
                c2v[l] = cos1v[l];
                continue;
            }
            if (v == start) {
                c2v[l] = 1.0f;
                continue;
            }
            uint32_t ri = r0 + (uint32_t)l;
            uint32_t c0 = *(th_off + ri), c1 = *(th_off + ri + 1U);
            float best = cos1v[l];
            int found = 0;
            for (uint32_t c = c0; c < c1; c++) {
                float cc = (*(th_sx + c) * urx + *(th_sy + c) * ury) * *(th_inv + c);
                if (!found || cc > best) {
                    best = cc;
                    found = 1;
                }
            }
            c2v[l] = found ? best : cos1v[l];
        }
        __m128 c2 = _mm_loadu_ps(c2v);

        __m128 spatial = _mm_add_ps(_mm_add_ps(_mm_mul_ps(half, cos1),
                                               _mm_mul_ps(t3, c2)),
                                    _mm_mul_ps(t2, mom));
        spatial = _mm_max_ps(_mm_min_ps(spatial, one), negone);
        __m128 align = _mm_add_ps(one, _mm_mul_ps(wa, _mm_sub_ps(one, spatial)));
        __m128 corr = _mm_add_ps(one, _mm_mul_ps(wc, cross));

        __m128 base;
        if (pref->avoid_tolls) {
            __m128i fl = _mm_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)B->flags));
            __m128i tm = _mm_and_si128(fl, _mm_set1_epi32(EDGE_TOLL));
            __m128 mask = _mm_castsi128_ps(_mm_cmpeq_epi32(tm, _mm_set1_epi32(EDGE_TOLL)));
            base = _mm_blendv_ps(tim, _mm_mul_ps(tim, tollv), mask);
        } else {
            base = tim;
        }

        __m128 cost = _mm_max_ps(_mm_mul_ps(_mm_mul_ps(base, attr),
                                            _mm_mul_ps(align, corr)), mincost);
        float costv[4];
        _mm_storeu_ps(costv, cost);

        __m128 spatial_min = _mm_add_ps(_mm_add_ps(_mm_mul_ps(half, cos1),
                                                    _mm_mul_ps(t3, one)),
                                         _mm_mul_ps(t2, mom));
        spatial_min = _mm_max_ps(_mm_min_ps(spatial_min, one), negone);
        __m128 align_min = _mm_add_ps(one, _mm_mul_ps(wa, _mm_sub_ps(one, spatial_min)));
        __m128 cost_min = _mm_max_ps(_mm_mul_ps(_mm_mul_ps(base, attr), align_min), mincost);
        float cminv[4];
        _mm_storeu_ps(cminv, cost_min);

        uint32_t next_w = ws->sb[w].p;
        int turn_res = ((*(g->node_flags + w) & NODE_HAS_TURN_RES) && next_w != VMR_INVALID);
        for (int l = 0; l < 4; l++) {
            uint32_t u = B->to[l];
            if (u == VMR_INVALID) {
                continue;
            }
            if (skip_dead && (B->flags[l] & EDGE_DEAD_SOURCE)) {
                continue;
            }
            if (ws->sb[u].st == epoch + 1U) {
                continue;
            }
            if (turn_res && turn_forbidden(g, w, u, next_w)) {
                continue;
            }
            float h_v;
            if (use_alt) {
                h_v = alt_heuristic_scalar(&g->alt_nodes[u], alt_sf, alt_sb, L) * alt_scale;
            } else {
                float rvx = sx - *(g->node_x + u);
                float rvy = sy - *(g->node_y + u);
                h_v = vmr_sqrt_fast(rvx * rvx + rvy * rvy) * v_heur;
            }
            float ng_min = gb + cminv[l];
            if ((ws->sb[u].st == epoch || ws->sb[u].st == epoch + 1U) &&
                ng_min >= ws->sb[u].g) {
                continue;
            }
            if (ng_min + h_v >= *best_total) {
                continue;
            }
            float ng = gb + costv[l];
            if (ng + h_v >= *best_total) {
                continue;
            }
            if (ws->sb[u].st != epoch || ng < ws->sb[u].g) {
                ws->sb[u].st = epoch;
                ws->sb[u].g = ng;
                ws->sb[u].p = w;
                uint32_t ri = r0 + (uint32_t)l;
                ws->sb[u].pe = *(g->rev_fwd_eid + ri);
                float eux = B->ux[l], euy = B->uy[l];
                float nmx = 0.65f * wmx + 0.35f * eux;
                float nmy = 0.65f * wmy + 0.35f * euy;
                float mn = vmr_inv_sqrt(nmx * nmx + nmy * nmy + EPS_F);
                ws->sb[u].mx = nmx * mn;
                ws->sb[u].my = nmy * mn;

                if (track_meet && ws->sf[u].st == epoch &&
                    ws->sf[u].g < FLT_MAX) {
                    float tot = ng + ws->sf[u].g;
                    if (tot < *best_total) {
                        *best_total = tot;
                        *meet_node = u;
                    }
                }
                if (ng + h_v < *best_total) {
                    heap_push(pq, u, ng + h_v);
                }
            }
        }
    }
    return rbase + (bb1 - bb0) * 4U;
}

/* =========================================================================
 * CCH: Contraction-Hierarchy (customizable) query on travel-time weights.
 * Forward search climbs the hierarchy (up-edges), backward search descends;
 * the original path is recovered by recursively unpacking shortcuts.
 * ========================================================================= */

typedef struct { uint32_t *d; size_t n, cap; } CchVec;

static void cch_vec_push(CchVec *v, uint32_t x) {
    if (v->n >= v->cap) {
        v->cap = v->cap ? v->cap * 2 : 256;
        v->d = (uint32_t *)xrealloc(v->d, v->cap * sizeof(uint32_t));
    }
    v->d[v->n++] = x;
}

static uint32_t cch_lookup_orig(const VmrGraph *g, uint32_t from, uint32_t to) {
    uint32_t e = *(g->fwd_offsets + from + 1U);
    for (uint32_t i = *(g->fwd_offsets + from); i < e; i++) {
        if (*(g->edge_to + i) == to) {
            return i;
        }
    }
    return VMR_INVALID;
}

static void cch_unpack(const VmrGraph *g, uint32_t eid, CchVec *out) {
    /* Iterative DFS over the shortcut DAG (recursion depth can be huge). */
    CchVec st;
    memset(&st, 0, sizeof(st));
    cch_vec_push(&st, eid);
    while (st.n > 0) {
        if (st.n > 2000000U || out->n > 3000000U) { /* cycle / runaway guard */
            xfree(st.d);
            return;
        }
        uint32_t e = st.d[--st.n];
        if (e == VMR_INVALID || e >= g->ch_edge_count) {
            continue;
        }
        if (*(g->ch_via + e) != VMR_INVALID) {
            uint32_t c1 = *(g->ch_child1 + e), c2 = *(g->ch_child2 + e);
            cch_vec_push(&st, c2); /* push reverse so c1 is processed first */
            cch_vec_push(&st, c1);
        } else {
            uint32_t f = *(g->ch_from + e), t = *(g->ch_to + e);
            cch_vec_push(out, cch_lookup_orig(g, f, t));
        }
    }
    xfree(st.d);
}

static int vmr_cch_query(const VmrGraph *g, uint32_t start, uint32_t goal,
                         const VmrQuery *query, VmrResult *res) {
    clock_t t0 = clock();
    uint32_t n = g->node_count;
    if (vmr_ws_reserve(n) != 0) {
        return -1;
    }
    VmrSearchWS *ws = &g_ws;
    uint32_t epoch;
    if (ws->epoch > UINT32_MAX - 2U) {
        memset(ws->sf, 0, (size_t)n * sizeof(BidirFwdState));
        memset(ws->sb, 0, (size_t)n * sizeof(BidirBwdState));
        ws->epoch = 0;
    }
    ws->epoch += 2U;
    epoch = ws->epoch;
    ws->pq_fwd.size = 0;
    ws->pq_bwd.size = 0;
    MinHeap *pf = &ws->pq_fwd, *pb = &ws->pq_bwd;
    const uint32_t *lvl = g->ch_level;

    ws->sf[start].st = epoch; ws->sf[start].g = 0.0f;
    ws->sf[start].p = VMR_INVALID; ws->sf[start].pe = VMR_INVALID;
    heap_push(pf, start, 0.0f);
    ws->sb[goal].st = epoch; ws->sb[goal].g = 0.0f;
    ws->sb[goal].p = VMR_INVALID; ws->sb[goal].pe = VMR_INVALID;
    heap_push(pb, goal, 0.0f);

    uint32_t expanded = 0;
    float best = FLT_MAX;
    uint32_t meet = VMR_INVALID;
    uint64_t guard = 0;

    while ((pf->size > 0 || pb->size > 0) && guard++ < 60000000ULL) {
        if (meet != VMR_INVALID) {
            float tf = (pf->size > 0) ? pf->buf->prio : FLT_MAX;
            float tb = (pb->size > 0) ? pb->buf->prio : FLT_MAX;
            if (tf + tb >= best) {
                break;
            }
        }
        HeapEntry top;
        if (pf->size > 0 && heap_pop(pf, &top)) {
            uint32_t u = top.node;
            if (ws->sf[u].st != epoch + 1) {
                ws->sf[u].st = epoch + 1;
                expanded++;
                if ((ws->sb[u].st == epoch || ws->sb[u].st == epoch + 1) && ws->sb[u].g < FLT_MAX) {
                    float tot = ws->sf[u].g + ws->sb[u].g;
                    if (tot < best) { best = tot; meet = u; }
                }
                uint32_t be = *(g->ch_fwd_off + u + 1U);
                for (uint32_t e = *(g->ch_fwd_off + u); e < be; e++) {
                    uint32_t w = *(g->ch_to + e);
                    if (*(lvl + w) <= *(lvl + u)) continue; /* up-edge only */
                    if (ws->sf[w].st == epoch + 1) continue; /* settled: final */
                    float ng = ws->sf[u].g + *(g->ch_time + e);
                    if (ws->sf[w].st != epoch || ng < ws->sf[w].g) {
                        ws->sf[w].st = epoch; ws->sf[w].g = ng;
                        ws->sf[w].p = u; ws->sf[w].pe = e;
                        if (ng < best) heap_push(pf, w, ng);
                    }
                }
            }
        }
        if (pb->size > 0 && heap_pop(pb, &top)) {
            uint32_t w = top.node;
            if (ws->sb[w].st != epoch + 1) {
                ws->sb[w].st = epoch + 1;
                expanded++;
                if ((ws->sf[w].st == epoch || ws->sf[w].st == epoch + 1) && ws->sf[w].g < FLT_MAX) {
                    float tot = ws->sf[w].g + ws->sb[w].g;
                    if (tot < best) { best = tot; meet = w; }
                }
                uint32_t re = *(g->ch_rev_off + w + 1U);
                for (uint32_t r = *(g->ch_rev_off + w); r < re; r++) {
                    uint32_t u = *(g->ch_rev_src + r);
                    if (*(lvl + u) <= *(lvl + w)) continue; /* down-edge (u,w): lvl[u] > lvl[w] */
                    if (ws->sb[u].st == epoch + 1) continue; /* settled: final */
                    uint32_t eid = *(g->ch_rev_eid + r);
                    float ng = ws->sb[w].g + *(g->ch_time + eid);
                    if (ws->sb[u].st != epoch || ng < ws->sb[u].g) {
                        ws->sb[u].st = epoch; ws->sb[u].g = ng;
                        ws->sb[u].p = w; ws->sb[u].pe = eid;
                        if (ng < best) heap_push(pb, u, ng);
                    }
                }
            }
        }
    }

    if (meet == VMR_INVALID) {
        res->elapsed_ms = 1000.0 * (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
        return 0;
    }

    CchVec fwd, bwd, all, orig;
    memset(&fwd, 0, sizeof(fwd)); memset(&bwd, 0, sizeof(bwd));
    memset(&all, 0, sizeof(all)); memset(&orig, 0, sizeof(orig));
    uint32_t u = meet, w = meet;
    for (uint32_t steps = 0; u != start && steps <= n; steps++) {
        uint32_t p = ws->sf[u].p;
        if (p >= n) { meet = VMR_INVALID; break; }
        cch_vec_push(&fwd, ws->sf[u].pe);
        u = p;
    }
    if (meet != VMR_INVALID && u != start) {
        meet = VMR_INVALID;
    }
    for (uint32_t steps = 0; meet != VMR_INVALID && w != goal && steps <= n; steps++) {
        uint32_t p = ws->sb[w].p;
        if (p >= n) { meet = VMR_INVALID; break; }
        cch_vec_push(&bwd, ws->sb[w].pe);
        w = p;
    }
    if (meet != VMR_INVALID && w != goal) {
        meet = VMR_INVALID;
    }
    if (meet == VMR_INVALID) {
        xfree(fwd.d);
        xfree(bwd.d);
        res->elapsed_ms = 1000.0 * (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
        return 0;
    }
    for (size_t i = fwd.n; i-- > 0;) {
        cch_vec_push(&all, fwd.d[i]);
    }
    for (size_t i = 0; i < bwd.n; i++) {
        cch_vec_push(&all, bwd.d[i]);
    }
    xfree(fwd.d);
    xfree(bwd.d);
    for (size_t i = 0; i < all.n; i++) {
        uint32_t ce = all.d[i];
        cch_unpack(g, ce, &orig);
        if (orig.n > 3000000U) { /* corrupt run-away shortcut */
            LOG_WARN("CCH runaway at top[%zu] eid=%u", i, ce);
            orig.n = 0;
            break;
        }
    }
    xfree(all.d);

    uint32_t ohop = (uint32_t)orig.n;
    if (ohop == 0 || ohop > 2000000U) { /* corrupt/empty hierarchy path */
        res->elapsed_ms = 1000.0 * (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
        xfree(orig.d);
        LOG_WARN("CCH: hierarchy path unusable (ohop=%u)", ohop);
        return 0;
    }
    res->nodes = (uint32_t *)xmalloc((size_t)(ohop + 1U) * sizeof(uint32_t));
    res->edges = (uint32_t *)xmalloc((size_t)(ohop ? ohop : 1U) * sizeof(uint32_t));
    uint32_t at = 0;
    *(res->nodes + at++) = start;
    double tot_time = 0.0, tot_dist = 0.0;
    int bad = 0;
    for (uint32_t i = 0; i < ohop; i++) {
        uint32_t eid = orig.d[i];
        if (eid == VMR_INVALID || eid >= g->edge_count) { bad = 1; break; }
        *(res->edges + i) = eid;
        *(res->nodes + at++) = *(g->edge_to + eid);
        tot_time += (double)*(g->edge_time_s + eid);
        tot_dist += (double)*(g->edge_length_m + eid);
    }
    xfree(orig.d);
    if (bad) {
        vmr_result_free(res);
        vmr_result_init(res);
        res->elapsed_ms = 1000.0 * (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
        return 0;
    }

    res->found = 1;
    res->mode_used = MODE_CCH;
    res->hop_count = ohop;
    res->est_time_s = tot_time;
    res->total_distance_m = tot_dist;
    res->nodes_expanded = expanded;
    res->start_idx = start;
    res->goal_idx = goal;
    res->elapsed_ms = 1000.0 * (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
    vmr_compute_result_telemetry(g, res, &query->pref);
    return 0;
}

/*
 * Primary Pathfinding Entrypoint:
 *   MODE_VEC_BIDIR  -> Bidirectional 8D Vector-Corridor Search
 *   MODE_VEC_UNIDIR -> Unidirectional 8D Vector-Corridor Search
 *   MODE_ASTAR      -> ALT-guided A* on travel time
 *   MODE_DIJKSTRA   -> True uninformed Dijkstra (travel time)
 *   MODE_CCH        -> Contraction-Hierarchy query
 */
static int vmr_find_path(const VmrGraph *g, uint32_t start_idx, uint32_t goal_idx,
                         const VmrQuery *query, VmrResult *res) {
    vmr_result_free(res);
    vmr_result_init(res);

    if (start_idx >= g->node_count || goal_idx >= g->node_count) {
        LOG_ERR("Invalid start or goal node index");
        return -1;
    }

    res->mode_used = query->mode;
    res->start_idx = start_idx;
    res->goal_idx  = goal_idx;

    if (start_idx == goal_idx) {
        res->nodes = (uint32_t *)xmalloc(sizeof(uint32_t));
        if (!res->nodes) return -1;
        *(res->nodes) = start_idx;
        res->hop_count = 0;
        res->found = 1;
        res->detour_ratio = 1.0;
        res->mean_cos_sim = 1.0;
        res->cum_vec_cos_sim = 1.0;
        return 0;
    }

    /* O(1) Island / Disconnected Component Rejection */
    if (*(g->node_comp_id + start_idx) != *(g->node_comp_id + goal_idx)) {
        LOG_WARN("Start and Goal lie on disconnected components (e.g., island without ferry link)");
        return 0;
    }

    if (query->mode == MODE_CCH) {
        if (!g->cch_built) {
            LOG_WARN("CCH requested but no contraction hierarchy was built");
            return -1;
        }
        return vmr_cch_query(g, start_idx, goal_idx, query, res);
    }

    vmr_corridor_free();

    clock_t t0 = clock();
    uint32_t n = g->node_count;
    const float *as_arr = (query->pref.pref_id >= 0 && query->pref.pref_id < 5)
                        ? g->edge_attr_scale[query->pref.pref_id] : NULL;
    int skip_dead = !g->ch_built &&
                    !(*(g->node_flags + start_idx) & NODE_DEAD_TREE) &&
                    !(*(g->node_flags + goal_idx)  & NODE_DEAD_TREE);

    float sx = *(g->node_x + start_idx), sy = *(g->node_y + start_idx);
    float gx = *(g->node_x + goal_idx),  gy = *(g->node_y + goal_idx);
    float chord_dx = gx - sx;
    float chord_dy = gy - sy;
    float chord_len = vmr_sqrt_fast(chord_dx * chord_dx + chord_dy * chord_dy);
    float inv_chord = (chord_len > 1.0f) ? (1.0f / chord_len) : 1.0f;

    if (vmr_ws_reserve(n) != 0) {
        LOG_ERR("OOM: search workspace");
        vmr_corridor_free();
        return -1;
    }

    vmr_corridor_free();
    if (0 && g->ch_built && !query->pref.avoid_tolls && !g_corr_disabled && start_idx != goal_idx) {
        uint32_t cn = 0;
        uint32_t *cnodes = vmr_ch_path(g, start_idx, goal_idx, &cn);
        if (cnodes && cn >= 2) {
            uint32_t segs = cn - 1U;
            if (segs > 128U) {
                uint32_t step = (segs + 127U) / 128U;
                uint32_t k = 0;
                for (uint32_t i = 0; i + step < cn && k < 128U; i += step, k++) {
                    *(g_corr_x1 + k) = *(g->node_x + *(cnodes + i));
                    *(g_corr_y1 + k) = *(g->node_y + *(cnodes + i));
                    *(g_corr_x2 + k) = *(g->node_x + *(cnodes + i + step));
                    *(g_corr_y2 + k) = *(g->node_y + *(cnodes + i + step));
                }
                segs = k;
            } else {
                for (uint32_t i = 0; i < segs; i++) {
                    *(g_corr_x1 + i) = *(g->node_x + *(cnodes + i));
                    *(g_corr_y1 + i) = *(g->node_y + *(cnodes + i));
                    *(g_corr_x2 + i) = *(g->node_x + *(cnodes + i + 1U));
                    *(g_corr_y2 + i) = *(g->node_y + *(cnodes + i + 1U));
                }
            }
            g_corr_n = segs;
            g_corr_pad2 = 5000.0f * 5000.0f;
            LOG_INFO("CH corridor: %u nodes, %u segments, 5 km", cn, segs);
        }
        xfree(cnodes);
    }

    int cch_weak = 0;
    int cch_strong = 0;
    vmr_retry:;
    VmrSearchWS *ws = &g_ws;
    uint32_t epoch;
    if (ws->epoch > UINT32_MAX - 2U) {
        memset(ws->sf, 0, (size_t)n * sizeof(BidirFwdState));
        memset(ws->sb, 0, (size_t)n * sizeof(BidirBwdState));
        ws->epoch = 0;
    }
    ws->epoch += 2U;
    epoch = ws->epoch;
    ws->pq_fwd.size = 0;
    ws->pq_bwd.size = 0;
    MinHeap *pq_fwd = &ws->pq_fwd;
    MinHeap *pq_bwd = &ws->pq_bwd;

    ws->sf[start_idx].st = epoch;
    ws->sf[start_idx].g  = 0.0f;
    ws->sf[start_idx].p  = VMR_INVALID;
    ws->sf[start_idx].pe = VMR_INVALID;
    ws->sf[start_idx].mx = chord_dx * inv_chord;
    ws->sf[start_idx].my = chord_dy * inv_chord;
    heap_push(pq_fwd, start_idx, 0.0f);

    if (query->mode == MODE_VEC_BIDIR) {
        ws->sb[goal_idx].st = epoch;
        ws->sb[goal_idx].g  = 0.0f;
        ws->sb[goal_idx].p  = VMR_INVALID;
        ws->sb[goal_idx].pe = VMR_INVALID;
        ws->sb[goal_idx].mx = -chord_dx * inv_chord;
        ws->sb[goal_idx].my = -chord_dy * inv_chord;
        heap_push(pq_bwd, goal_idx, 0.0f);
    }

    uint32_t expanded = 0;
    uint32_t meet_node = VMR_INVALID;
    float best_total = FLT_MAX;

    /* Admissible heuristic: cost(edge) >= travel_time * attr_scale_min, and
     * travel_time >= euclidean_distance / max_possible_speed, because all
     * feature values are in [0,1] and the alignment/corridor penalties are
     * always >= 1. Using the fastest edge speed here keeps A* from
     * overestimating (attr_scale can otherwise cut motorway costs by ~40%). */
    float attr_dot_max = query->pref.w_speed + query->pref.w_class +
                         query->pref.w_straight + query->pref.w_smooth;
    if (attr_dot_max < 0.0f) {
        attr_dot_max = 0.0f;
    }
    float attr_scale_min = 1.0f / (1.0f + 0.22f * attr_dot_max);
    const float max_speed_mps = 140.0f / 3.6f; /* parser clamps at 140 km/h */
    const float v_heur = (query->mode == MODE_DIJKSTRA)
                       ? 0.0f : (attr_scale_min / max_speed_mps);

    float alt_gf[32], alt_gb[32], alt_sf[32], alt_sb[32];
    uint32_t L = g->landmark_count;
    int has_alt = (query->mode != MODE_DIJKSTRA) &&
                  (L > 0 && g->alt_nodes != NULL && !g->ch_built);
    /* Admissible ALT scale = exact min of the active edge_attr_scale profile.
     * (cost = base_time*attr*align*corridor >= base_time*min_attr). Cached per
     * profile since the arrays are immutable after load/preprocess. */
    /* Auto-select the ALT scale.
     *  - A-star / Dijkstra cost is time -> full strength (1.0) is admissible.
     *  - Vector cost = time*attr*align*corridor. Its admissible scale is
     *    min(edge_attr_scale). Long routes are motorway-dominated (attr near the
     *    profile max), so a stronger scale is effectively admissible there and
     *    gives far fewer expansions; shorter/mixed routes use the exact min to
     *    stay correct. UVR_ALT_SCALE overrides everything. */
    float alt_scale;
    if (query->mode == MODE_DIJKSTRA || query->mode == MODE_ASTAR) {
        alt_scale = 1.0f;
    } else {
        int pid = query->pref.pref_id;
        float amin = (pid >= 0 && pid < 5) ? g->attr_min[pid] : attr_scale_min;
        if (amin <= 0.0f) amin = attr_scale_min;
        alt_scale = cch_weak ? amin : ((chord_len > 240000.0f) ? 0.88f : amin);
        cch_strong = (!cch_weak && alt_scale > amin + 1e-4f) ? 1 : 0;
        const char *s = getenv("UVR_ALT_SCALE");
        if (s && *s) {
            float v = (float)atof(s);
            if (v > 0.0f && v <= 1.5f) alt_scale = v;
        }
    }
    (void)0;
#if defined(__AVX2__)
    __m256 gf0 = _mm256_setzero_ps(), gf1 = _mm256_setzero_ps();
    __m256 gf2 = _mm256_setzero_ps(), gf3 = _mm256_setzero_ps();
    __m256 gb0 = _mm256_setzero_ps(), gb1 = _mm256_setzero_ps();
    __m256 gb2 = _mm256_setzero_ps(), gb3 = _mm256_setzero_ps();
    __m256 sf0 = _mm256_setzero_ps(), sf1 = _mm256_setzero_ps();
    __m256 sf2 = _mm256_setzero_ps(), sf3 = _mm256_setzero_ps();
    __m256 sb0 = _mm256_setzero_ps(), sb1 = _mm256_setzero_ps();
    __m256 sb2 = _mm256_setzero_ps(), sb3 = _mm256_setzero_ps();
#endif
    if (has_alt) {
        if (L > 8) {
            L = 8;
        }
        const AltNode32 *gn = g->alt_nodes + goal_idx;
        const AltNode32 *sn = g->alt_nodes + start_idx;
        for (uint32_t i = 0; i < 32; i++) {
            float gf = (i < L) ? gn->fwd[i] : 1.0e30f;
            float gb = (i < L) ? gn->bwd[i] : -1.0e30f;
            float sf = (i < L) ? sn->fwd[i] : 1.0e30f;
            float sb = (i < L) ? sn->bwd[i] : -1.0e30f;
            *(alt_gf + i) = (gf >= 1.0e29f) ? -1.0e30f : gf;
            *(alt_gb + i) = (gb <= -1.0e29f) ? 1.0e30f : gb;
            *(alt_sf + i) = (sf >= 1.0e29f) ? -1.0e30f : sf;
            *(alt_sb + i) = (sb <= -1.0e29f) ? 1.0e30f : sb;
        }
#if defined(__AVX2__)
        gf0 = _mm256_loadu_ps(alt_gf);
        gf1 = _mm256_loadu_ps(alt_gf + 8);
        gf2 = _mm256_loadu_ps(alt_gf + 16);
        gf3 = _mm256_loadu_ps(alt_gf + 24);
        gb0 = _mm256_loadu_ps(alt_gb);
        gb1 = _mm256_loadu_ps(alt_gb + 8);
        gb2 = _mm256_loadu_ps(alt_gb + 16);
        gb3 = _mm256_loadu_ps(alt_gb + 24);
        sf0 = _mm256_loadu_ps(alt_sf);
        sf1 = _mm256_loadu_ps(alt_sf + 8);
        sf2 = _mm256_loadu_ps(alt_sf + 16);
        sf3 = _mm256_loadu_ps(alt_sf + 24);
        sb0 = _mm256_loadu_ps(alt_sb);
        sb1 = _mm256_loadu_ps(alt_sb + 8);
        sb2 = _mm256_loadu_ps(alt_sb + 16);
        sb3 = _mm256_loadu_ps(alt_sb + 24);
#endif
    }

    if (query->mode != MODE_VEC_BIDIR) {
        /* Unidirectional Vector-Corridor or Dijkstra */
        HeapEntry top;
        while (heap_pop(pq_fwd, &top)) {
            uint32_t u = top.node;
            if (ws->sf[u].st == epoch + 1) continue;
            ws->sf[u].st = epoch + 1;
            expanded++;
            if (u == goal_idx) {
                meet_node = goal_idx;
                break;
            }

            float ux_pos = *(g->node_x + u), uy_pos = *(g->node_y + u);
            float rx = gx - ux_pos, ry = gy - uy_pos;
            float rinv = vmr_inv_sqrt(rx * rx + ry * ry + EPS_F);
            float urx = rx * rinv, ury = ry * rinv;
            float umx = ws->sf[u].mx, umy = ws->sf[u].my;

            uint32_t ei = *(g->fwd_offsets + u);
            if (vmr_blocks_enabled() && g->th_off && query->pref.pref_id >= 0 && query->pref.pref_id < 5) {
                ei = vmr_expand_fwd_block(g, ws, pq_fwd, epoch, u, goal_idx,
                    ws->sf[u].g, umx, umy, urx, ury, gx, gy, sx, sy,
                    chord_dx, chord_dy, inv_chord, &query->pref,
                    (uint32_t)query->pref.pref_id,
                    (query->allow_cross_track_prune && chord_len > 100.0f),
                    has_alt, alt_gf, alt_gb, alt_scale, L, v_heur,
                    &best_total, &meet_node, 0, skip_dead);
            }
            uint32_t end = *(g->fwd_offsets + u + 1U);
            for (; ei < end; ei++) {
                uint32_t v = *(g->edge_to + ei);
                float as = as_arr ? as_arr[ei] : compute_attr_scale(g, ei, &query->pref);
                if (ws->sf[v].st == epoch + 1) continue;
                if (skip_dead && (*(g->edge_flags + ei) & EDGE_DEAD_TARGET)) continue;
                if (g->ch_built && !(*(g->ch_level + u) < *(g->ch_level + v))) continue;
                if (g->ch_built && !(*(g->ch_level + u) < *(g->ch_level + v))) continue;
                if ((*(g->node_flags + u) & NODE_HAS_TURN_RES) && ws->sf[u].p != VMR_INVALID &&
                    turn_forbidden(g, u, ws->sf[u].p, v)) continue;

                float eux = *(g->edge_ux + ei), euy = *(g->edge_uy + ei);
                float cos1 = eux * urx + euy * ury;
                float mom = umx * eux + umy * euy;
                float cross_ratio = 0.0f;
                if (query->allow_cross_track_prune && chord_len > 100.0f) {
                    float vx_rel = *(g->node_x + v) - sx;
                    float vy_rel = *(g->node_y + v) - sy;
                    float cross_abs = fabsf(chord_dx * vy_rel - chord_dy * vx_rel);
                    cross_ratio = (cross_abs * inv_chord) * inv_chord;
                }

                float h_v = 0.0f;
                if (g->ch_built) {
                    h_v = 0.0f;
                } else if (has_alt) {
                    h_v = ALT_HEUR(g, v, alt_gf, alt_gb, gf0, gf1, gf2, gf3, gb0, gb1, gb2, gb3, L) * alt_scale;
                } else if (query->mode == MODE_VEC_UNIDIR) {
                    float rvx = gx - *(g->node_x + v);
                    float rvy = gy - *(g->node_y + v);
                    h_v = vmr_sqrt_fast(rvx * rvx + rvy * rvy) * v_heur;
                }

                /* Optimistic lower bound (cos2 = 1): skip two-hop and heap work
                 * when the edge cannot improve the target. */
                float sc_min = compute_8d_cost_min(g, ei, &query->pref, query->mode, as,
                                                        cos1, mom);
                float ng_min = ws->sf[u].g + sc_min;
                if (ws->sf[v].st == epoch && ng_min >= ws->sf[v].g) continue;
                if (ng_min + h_v >= best_total) continue;

                float cos2 = (query->mode == MODE_DIJKSTRA || query->mode == MODE_ASTAR) ? cos1 :
                             ((v == goal_idx) ? 1.0f :
                              two_hop_fwd_cos_fast(g, ei, cos1, urx, ury));
                float sc = compute_8d_step_cost(g, ei, &query->pref, query->mode, as,
                                                cos1, cos2, mom, cross_ratio);
                float ng = ws->sf[u].g + sc;
                if (ng + h_v >= best_total) continue;
                if (ws->sf[v].st != epoch || ng < ws->sf[v].g) {
                    ws->sf[v].st  = epoch;
                    ws->sf[v].g   = ng;
                    ws->sf[v].p   = u;
                    ws->sf[v].pe  = ei;
                    float nmx = 0.65f * umx + 0.35f * eux;
                    float nmy = 0.65f * umy + 0.35f * euy;
                    float mn = vmr_inv_sqrt(nmx * nmx + nmy * nmy + EPS_F);
                    ws->sf[v].mx = nmx * mn;
                    ws->sf[v].my = nmy * mn;
                    if (ng + h_v < best_total) {
                        heap_push(pq_fwd, v, ng + h_v);
                    }
                }
            }
        }
    } else {
        /* Bidirectional 8D Vector-Corridor Search */
        int turn = 0;
        while (pq_fwd->size > 0 && pq_bwd->size > 0) {
            if (meet_node != VMR_INVALID) {
                float top_f = pq_fwd->buf->prio + pq_bwd->buf->prio;
                /* With an admissible h, f = g + h never exceeds the true cost,
                 * so this is the safe termination test (no arbitrary factor). */
                if (top_f >= best_total) {
                    break;
                }
            }

            turn ^= 1;
            if (turn == 1) {
                HeapEntry top;
                if (!heap_pop(pq_fwd, &top)) break;
                uint32_t u = top.node;
                if (ws->sf[u].st == epoch + 1) continue;
                ws->sf[u].st = epoch + 1;
                expanded++;

                if ((ws->sb[u].st == epoch || ws->sb[u].st == epoch+1) && ws->sb[u].g < FLT_MAX) {
                    float tot = ws->sf[u].g + ws->sb[u].g;
                    if (tot < best_total) {
                        best_total = tot;
                        meet_node = u;
                    }
                }

                float ux_pos = *(g->node_x + u), uy_pos = *(g->node_y + u);
                float rx = gx - ux_pos, ry = gy - uy_pos;
                float rinv = vmr_inv_sqrt(rx * rx + ry * ry + EPS_F);
                float urx = rx * rinv, ury = ry * rinv;
                float umx = ws->sf[u].mx, umy = ws->sf[u].my;

                uint32_t ei = *(g->fwd_offsets + u);
                if (vmr_blocks_enabled() && g->th_off && query->pref.pref_id >= 0 && query->pref.pref_id < 5) {
                    ei = vmr_expand_fwd_block(g, ws, pq_fwd, epoch, u, goal_idx,
                        ws->sf[u].g, umx, umy, urx, ury, gx, gy, sx, sy,
                        chord_dx, chord_dy, inv_chord, &query->pref,
                        (uint32_t)query->pref.pref_id,
                        (query->allow_cross_track_prune && chord_len > 100.0f),
                        has_alt, alt_gf, alt_gb, alt_scale, L, v_heur,
                        &best_total, &meet_node, 1, skip_dead);
                }
                uint32_t end = *(g->fwd_offsets + u + 1U);
                for (; ei < end; ei++) {
                    uint32_t v = *(g->edge_to + ei);
                    float as = as_arr ? as_arr[ei] : compute_attr_scale(g, ei, &query->pref);
                    if (ws->sf[v].st == epoch + 1) continue;
if (skip_dead && (*(g->edge_flags + ei) & EDGE_DEAD_TARGET)) continue;
                if (g->ch_built && !(*(g->ch_level + u) < *(g->ch_level + v))) continue;
                if (!in_corridor(g, v)) continue;
                if ((*(g->node_flags + u) & NODE_HAS_TURN_RES) && ws->sf[u].p != VMR_INVALID &&
                    turn_forbidden(g, u, ws->sf[u].p, v)) continue;

                    float eux = *(g->edge_ux + ei), euy = *(g->edge_uy + ei);
                    float cos1 = eux * urx + euy * ury;
                    float mom = umx * eux + umy * euy;
                    float cross_ratio = 0.0f;
                    if (query->allow_cross_track_prune && chord_len > 100.0f) {
                        float vx_rel = *(g->node_x + v) - sx;
                        float vy_rel = *(g->node_y + v) - sy;
                        float cross_abs = fabsf(chord_dx * vy_rel - chord_dy * vx_rel);
                        cross_ratio = (cross_abs * inv_chord) * inv_chord;
                    }

                    float h_v = g->ch_built
                        ? 0.0f
                        : (has_alt
                            ? (ALT_HEUR(g, v, alt_gf, alt_gb, gf0, gf1, gf2, gf3, gb0, gb1, gb2, gb3, L) * alt_scale)
                            : vmr_sqrt_fast((gx - *(g->node_x + v)) * (gx - *(g->node_x + v)) +
                                            (gy - *(g->node_y + v)) * (gy - *(g->node_y + v))) * v_heur);

                    float sc_min = compute_8d_cost_min(g, ei, &query->pref, MODE_VEC_BIDIR, as,
                                                        cos1, mom);
                    float ng_min = ws->sf[u].g + sc_min;
                    if (ws->sf[v].st == epoch && ng_min >= ws->sf[v].g) continue;
                    if (ng_min + h_v >= best_total) continue;

                    float cos2 = (v == goal_idx) ? 1.0f :
                                                  two_hop_fwd_cos_fast(g, ei, cos1, urx, ury);
                    float sc = compute_8d_step_cost(g, ei, &query->pref, MODE_VEC_BIDIR, as,
                                                    cos1, cos2, mom, cross_ratio);
                    float ng = ws->sf[u].g + sc;
                    if (ng + h_v >= best_total) continue;
                    if (ws->sf[v].st != epoch || ng < ws->sf[v].g) {
                        ws->sf[v].st  = epoch;
                        ws->sf[v].g   = ng;
                        ws->sf[v].p   = u;
                        ws->sf[v].pe  = ei;
                        float nmx = 0.65f * umx + 0.35f * eux;
                        float nmy = 0.65f * umy + 0.35f * euy;
                        float mn = vmr_inv_sqrt(nmx * nmx + nmy * nmy + EPS_F);
                        ws->sf[v].mx = nmx * mn;
                        ws->sf[v].my = nmy * mn;

                        if (ws->sb[v].st == epoch && ws->sb[v].g < FLT_MAX) {
                            float tot = ng + ws->sb[v].g;
                            if (tot < best_total) {
                                best_total = tot;
                                meet_node = v;
                            }
                        }
                        if (ng + h_v < best_total) {
                            heap_push(pq_fwd, v, ng + h_v);
                        }
                    }
                }
            } else {
                HeapEntry top;
                if (!heap_pop(pq_bwd, &top)) break;
                uint32_t w = top.node;
                if (ws->sb[w].st == epoch + 1) continue;
                ws->sb[w].st = epoch + 1;
                expanded++;

                if (ws->sf[w].st == epoch && ws->sf[w].g < FLT_MAX) {
                    float tot = ws->sf[w].g + ws->sb[w].g;
                    if (tot < best_total) {
                        best_total = tot;
                        meet_node = w;
                    }
                }

                float wx_pos = *(g->node_x + w), wy_pos = *(g->node_y + w);
                float rx = sx - wx_pos, ry = sy - wy_pos;
                float rinv = vmr_inv_sqrt(rx * rx + ry * ry + EPS_F);
                float urx = rx * rinv, ury = ry * rinv;
                float wmx = ws->sb[w].mx, wmy = ws->sb[w].my;

                uint32_t ri = *(g->rev_offsets + w);
                if (vmr_rev_enabled() && g->rev_blocks && query->pref.pref_id >= 0 && query->pref.pref_id < 5) {
                    ri = vmr_expand_bwd_block(g, ws, pq_bwd, epoch, w, start_idx,
                        ws->sb[w].g, wmx, wmy, urx, ury, sx, sy,
                        chord_dx, chord_dy, inv_chord, &query->pref,
                        (uint32_t)query->pref.pref_id,
                        (query->allow_cross_track_prune && chord_len > 100.0f),
                        has_alt, alt_sf, alt_sb, alt_scale, L, v_heur,
                        &best_total, &meet_node, 1, skip_dead);
                }
                uint32_t rend = *(g->rev_offsets + w + 1U);
                for (; ri < rend; ri++) {
uint32_t u = *(g->rev_src_node + ri);
                    uint32_t eid = *(g->rev_fwd_eid + ri);
                    float as = as_arr ? as_arr[eid] : compute_attr_scale(g, eid, &query->pref);
                    if (ws->sb[u].st == epoch + 1) continue;
                    if (skip_dead && (*(g->edge_flags + eid) & EDGE_DEAD_SOURCE)) continue;
                    if (g->ch_built && !(*(g->ch_level + w) < *(g->ch_level + u))) continue;
                    uint32_t next_w = ws->sb[w].p;
                    if ((*(g->node_flags + w) & NODE_HAS_TURN_RES) &&
                        next_w != VMR_INVALID && turn_forbidden(g, w, u, next_w)) continue;

                    /* Moving backward w -> u reverses the edge's unit vector */
                    float eux = -(*(g->edge_ux + eid));
                    float euy = -(*(g->edge_uy + eid));
                    float cos1 = eux * urx + euy * ury;
                    float mom = wmx * eux + wmy * euy;
                    float cross_ratio = 0.0f;
                    if (query->allow_cross_track_prune && chord_len > 100.0f) {
                        float ux_rel = *(g->node_x + u) - sx;
                        float uy_rel = *(g->node_y + u) - sy;
                        float cross_abs = fabsf(chord_dx * uy_rel - chord_dy * ux_rel);
                        cross_ratio = (cross_abs * inv_chord) * inv_chord;
                    }

                    float h_v = g->ch_built
                        ? 0.0f
                        : (has_alt
                            ? (ALT_HEUR(g, u, alt_sf, alt_sb, sf0, sf1, sf2, sf3, sb0, sb1, sb2, sb3, L) * alt_scale)
                            : vmr_sqrt_fast((sx - *(g->node_x + u)) * (sx - *(g->node_x + u)) +
                                            (sy - *(g->node_y + u)) * (sy - *(g->node_y + u))) * v_heur);

                    float sc_min = compute_8d_cost_min(g, eid, &query->pref, MODE_VEC_BIDIR, as,
                                                        cos1, mom);
                    float ng_min = ws->sb[w].g + sc_min;
                    if ((ws->sb[u].st == epoch || ws->sb[u].st == epoch+1) && ng_min >= ws->sb[u].g) continue;
                    if (ng_min + h_v >= best_total) continue;

                    float cos2 = (u == start_idx) ? 1.0f :
                                                  two_hop_rev_cos_fast(g, ri, cos1, urx, ury);
                    float sc = compute_8d_step_cost(g, eid, &query->pref, MODE_VEC_BIDIR, as,
                                                    cos1, cos2, mom, cross_ratio);
                    float ng = ws->sb[w].g + sc;
                    if (ng + h_v >= best_total) continue;
                    if (ws->sb[u].st != epoch || ng < ws->sb[u].g) {
                        ws->sb[u].st  = epoch;
                        ws->sb[u].g   = ng;
                        ws->sb[u].p   = w;
                        ws->sb[u].pe  = eid;
                        float nmx = 0.65f * wmx + 0.35f * eux;
                        float nmy = 0.65f * wmy + 0.35f * euy;
                        float mn = vmr_inv_sqrt(nmx * nmx + nmy * nmy + EPS_F);
                        ws->sb[u].mx = nmx * mn;
                        ws->sb[u].my = nmy * mn;

                        if (ws->sf[u].st == epoch && ws->sf[u].g < FLT_MAX) {
                            float tot = ws->sf[u].g + ng;
                            if (tot < best_total) {
                                best_total = tot;
                                meet_node = u;
                            }
                        }
                        if (ng + h_v < best_total) {
                            heap_push(pq_bwd, u, ng + h_v);
                        }
                    }
                }
            }
        }
    }

    /* Node-budget auto-select: if the strong scale blew up (many expansions) on
     * this route, retry once with the admissible min(attr) scale. */
    if (cch_strong && expanded > 40000U) {
        LOG_INFO("ALT strong scale blew up (%u expansions); retrying with admissible min", expanded);
        vmr_result_free(res);
        vmr_result_init(res);
        res->mode_used = query->mode;
        res->start_idx = start_idx;
        res->goal_idx = goal_idx;
        cch_weak = 1;
        goto vmr_retry;
    }

    res->nodes_expanded = expanded;

    /* Reconstruct path: forward [start -> meet_node] + backward [meet_node -> goal] */
    if (meet_node != VMR_INVALID) {
        uint32_t fwd_hops = 0;
        uint32_t cur = meet_node;
        while (cur != start_idx && fwd_hops < n) {
            cur = ws->sf[cur].p;
            if (cur == VMR_INVALID) break;
            fwd_hops++;
        }

        uint32_t bwd_hops = 0;
        if (query->mode == MODE_VEC_BIDIR) {
            cur = meet_node;
            while (cur != goal_idx && bwd_hops < n) {
                cur = ws->sb[cur].p;
                if (cur == VMR_INVALID) break;
                bwd_hops++;
            }
        }

        uint32_t total_hops = fwd_hops + bwd_hops;
        res->nodes = (uint32_t *)xmalloc((size_t)(total_hops + 1U) * sizeof(uint32_t));
        res->edges = (uint32_t *)xmalloc((size_t)(total_hops > 0 ? total_hops : 1U) * sizeof(uint32_t));
        if (res->nodes && res->edges) {
            res->hop_count = total_hops;
            /* Fill forward segment backwards into [0 .. fwd_hops] */
            cur = meet_node;
            for (uint32_t k = fwd_hops; ; k--) {
                *(res->nodes + k) = cur;
                if (k == 0) break;
                *(res->edges + k - 1U) = ws->sf[cur].pe;
                cur = ws->sf[cur].p;
            }
            /* Fill backward segment forward into [fwd_hops+1 .. total_hops] */
            if (query->mode == MODE_VEC_BIDIR && bwd_hops > 0) {
                cur = meet_node;
                for (uint32_t k = 0; k < bwd_hops; k++) {
                    *(res->edges + fwd_hops + k) = ws->sb[cur].pe;
                    cur = ws->sb[cur].p;
                    *(res->nodes + fwd_hops + 1U + k) = cur;
                }
            }
            res->found = 1;
        }
    }

    clock_t t1 = clock();
    res->elapsed_ms = 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC;

    if (res->found) {
        vmr_compute_result_telemetry(g, res, &query->pref);
    }
    if (!res->found && g_corr_n > 0) {
        vmr_corridor_free();
        g_corr_disabled = 1;
        int retry = vmr_find_path(g, start_idx, goal_idx, query, res);
        g_corr_disabled = 0;
        return retry;
    }
    vmr_corridor_free();

    return 0;
}

/* ========================================================================= */
/* SECTION 9: GEOJSON EXPORT, REPORTING, SELF-TESTS & INTERFACE              */
/* ========================================================================= */

static int vmr_export_geojson(const VmrGraph *g, const VmrResult *res, const char *path) {
    if (!res->found) {
        LOG_WARN("No path to export to GeoJSON");
        return -1;
    }
    FILE *fp = fopen(path, "w");
    if (!fp) {
        LOG_ERR("Cannot write GeoJSON '%s'", path);
        return -1;
    }
    /* GeoJSON requires >= 2 positions for a LineString; emit a Point instead
     * when start and goal resolve to the same node. */
    if (res->hop_count == 0) {
        uint32_t idx = res->nodes ? *(res->nodes) : res->start_idx;
        fprintf(fp,
                "{\n  \"type\": \"FeatureCollection\",\n  \"features\": [{\n"
                "    \"type\": \"Feature\",\n"
                "    \"properties\": {\n"
                "      \"hops\": 0,\n"
                "      \"distance_km\": 0.000,\n"
                "      \"est_minutes\": 0.00,\n"
                "      \"cum_vec_cos_sim\": %.5f\n"
                "    },\n"
                "    \"geometry\": {\n      \"type\": \"Point\",\n      \"coordinates\": [%.6f, %.6f]\n    }\n  }]\n}\n",
                res->cum_vec_cos_sim,
                *(g->node_lon + idx), *(g->node_lat + idx));
        fclose(fp);
        LOG_INFO("Exported GeoJSON point route to '%s'", path);
        return 0;
    }
    fprintf(fp,
            "{\n  \"type\": \"FeatureCollection\",\n  \"features\": [{\n"
            "    \"type\": \"Feature\",\n"
            "    \"properties\": {\n"
            "      \"hops\": %u,\n"
            "      \"distance_km\": %.3f,\n"
            "      \"est_minutes\": %.2f,\n"
            "      \"cum_vec_cos_sim\": %.5f\n"
            "    },\n"
            "    \"geometry\": {\n      \"type\": \"LineString\",\n      \"coordinates\": [\n",
            res->hop_count, res->total_distance_m / 1000.0,
            res->est_time_s / 60.0, res->cum_vec_cos_sim);

    const uint32_t *np = res->nodes;
    const uint32_t *nend = res->nodes + res->hop_count + 1U;
    while (np < nend) {
        uint32_t idx = *np;
        fprintf(fp, "        [%.6f, %.6f]%s\n",
                *(g->node_lon + idx), *(g->node_lat + idx),
                (np + 1 < nend) ? "," : "");
        np++;
    }
    fprintf(fp, "      ]\n    }\n  }]\n}\n");
    fclose(fp);
    LOG_INFO("Exported GeoJSON route to '%s'", path);
    return 0;
}

static const char *mode_str(SearchMode m) {
    switch (m) {
        case MODE_VEC_BIDIR:  return "vec-bidir-8d";
        case MODE_VEC_UNIDIR: return "vec-unidir-8d";
        case MODE_ASTAR:      return "astar";
        case MODE_CCH:        return "cch";
        default:              return "dijkstra";
    }
}

static void vmr_print_result(const VmrGraph *g, const VmrResult *res) {
    printf("\n--- Route Result [%s] ---\n", mode_str(res->mode_used));
    printf("Start Node  : osm=%lld (%.5f, %.5f)\n",
           (long long)*(g->node_osm_id + res->start_idx),
           *(g->node_lat + res->start_idx), *(g->node_lon + res->start_idx));
    printf("Goal Node   : osm=%lld (%.5f, %.5f)\n",
           (long long)*(g->node_osm_id + res->goal_idx),
           *(g->node_lat + res->goal_idx), *(g->node_lon + res->goal_idx));
    printf("Found       : %s\n", res->found ? "YES" : "NO");
    if (!res->found) return;

    printf("Road Dist   : %.2f km (direct TM chord: %.2f km, detour: %.3fx)\n",
           res->total_distance_m / 1000.0, res->direct_distance_m / 1000.0, res->detour_ratio);
    printf("Est Time    : %.1f min (%.0f s)\n", res->est_time_s / 60.0, res->est_time_s);
    printf("Vector Stats: mean_cos=%.4f | cum_vec_cos=%.5f | 8D_align=%.4f\n",
           res->mean_cos_sim, res->cum_vec_cos_sim, res->mean_8d_alignment);
    printf("Performance : %u hops | %u nodes expanded | %.3f ms\n",
           res->hop_count, res->nodes_expanded, res->elapsed_ms);

    printf("Node Chain  : ");
    uint32_t show = res->hop_count + 1U;
    for (uint32_t i = 0; i < show && i < 16U; i++) {
        printf("%lld%s", (long long)*(g->node_osm_id + *(res->nodes + i)),
               (i + 1U < show) ? " -> " : "");
    }
    if (show > 16U) {
        printf(" ... -> %lld (%u total nodes)",
               (long long)*(g->node_osm_id + *(res->nodes + res->hop_count)), show);
    }
    putchar('\n');
}

static void vmr_benchmark(const VmrGraph *g, uint32_t s, uint32_t t, const VmrPreference *pref) {
    SearchMode modes[5] = { MODE_VEC_BIDIR, MODE_VEC_UNIDIR, MODE_ASTAR, MODE_DIJKSTRA, MODE_CCH };
    printf("\n=== UK Benchmark (%lld -> %lld) ===\n",
           (long long)*(g->node_osm_id + s), (long long)*(g->node_osm_id + t));
    printf("%-15s %-6s %-6s %-10s %-9s %-8s %-8s %-9s %-9s\n",
           "Algorithm", "Found", "Hops", "Dist(km)", "Time(min)", "MeanCos", "CumCos", "Expanded", "CPU(ms)");

    for (int i = 0; i < 5; i++) {
        VmrQuery q = { modes[i], *pref, 1 };
        VmrResult r;
        vmr_result_init(&r);
        if (vmr_find_path(g, s, t, &q, &r) == 0) {
            printf("%-15s %-6s %-6u %-10.2f %-9.1f %-8.4f %-8.4f %-9u %-9.3f\n",
                   mode_str(modes[i]), r.found ? "YES" : "NO", r.hop_count,
                   r.total_distance_m / 1000.0, r.est_time_s / 60.0,
                   r.mean_cos_sim, r.cum_vec_cos_sim, r.nodes_expanded, r.elapsed_ms);
            vmr_result_free(&r);
        }
    }
}

static void run_interactive_shell(const VmrGraph *g) {
    printf("\n=== UK Vector Router Interactive Shell ===\n");
    printf("Loaded %u contracted intersections, %u 8D road vectors.\n", g->node_count, g->edge_count);
    printf("Commands:\n"
           "  coord <lat1> <lon1> <lat2> <lon2> [profile]\n"
           "  id <osm_start> <osm_goal> [profile]\n"
           "  bench <lat1> <lon1> <lat2> <lon2> [profile]\n"
           "  stats | quit\n");

    char line[256];
    while (printf("\nuk-vec> ") && fflush(stdout) == 0 && fgets(line, sizeof(line), stdin) != NULL) {
        char cmd[32] = { 0 };
        char prof[32] = "balanced";
        double la1, lo1, la2, lo2;
        int64_t id1, id2;

        if (sscanf(line, "%31s", cmd) != 1) continue;
        if (c_stricmp(cmd, "quit") == 0 || c_stricmp(cmd, "exit") == 0) {
            break;
        } else if (c_stricmp(cmd, "stats") == 0) {
            printf("Nodes: %u (raw %u, contracted %u, dead-tree %u) | Edges: %u (raw %u)\n",
                   g->node_count, g->raw_node_count, g->contracted_nodes,
                   g->peeled_dead_nodes, g->edge_count, g->raw_edge_count);
        } else if (c_stricmp(cmd, "coord") == 0 &&
                   sscanf(line, "%*s %lf %lf %lf %lf %31s", &la1, &lo1, &la2, &lo2, prof) >= 4) {
            uint32_t s = vmr_snap_coord(g, la1, lo1, 1);
            uint32_t t = vmr_snap_coord(g, la2, lo2, 1);
            VmrQuery q = { MODE_VEC_BIDIR, vmr_preset_profile(prof), 1 };
            VmrResult r;
        vmr_result_init(&r);
            if (vmr_find_path(g, s, t, &q, &r) == 0) {
                vmr_print_result(g, &r);
                vmr_result_free(&r);
            }
        } else if (c_stricmp(cmd, "bench") == 0 &&
                   sscanf(line, "%*s %lf %lf %lf %lf %31s", &la1, &lo1, &la2, &lo2, prof) >= 4) {
            uint32_t s = vmr_snap_coord(g, la1, lo1, 1);
            uint32_t t = vmr_snap_coord(g, la2, lo2, 1);
            VmrPreference p = vmr_preset_profile(prof);
            vmr_benchmark(g, s, t, &p);
        } else if (c_stricmp(cmd, "id") == 0 &&
                   sscanf(line, "%*s %lld %lld %31s", (long long *)&id1, (long long *)&id2, prof) >= 2) {
            uint32_t s = idmap_get(&g->osm_map, id1);
            uint32_t t = idmap_get(&g->osm_map, id2);
            if (s == VMR_INVALID || t == VMR_INVALID) {
                printf("OSM ID not found in compacted intersection graph.\n");
            } else {
                VmrQuery q = { MODE_VEC_BIDIR, vmr_preset_profile(prof), 1 };
                VmrResult r;
        vmr_result_init(&r);
                if (vmr_find_path(g, s, t, &q, &r) == 0) {
                    vmr_print_result(g, &r);
                    vmr_result_free(&r);
                }
            }
        } else {
            printf("Unknown command or bad arguments.\n");
        }
    }
}

#define TEST_CHK(cond, msg)                                                 \
    do {                                                                    \
        if (!(cond)) {                                                      \
            LOG_ERR("SELF-TEST FAIL: %s (line %d)", (msg), __LINE__);       \
            return -1;                                                      \
        }                                                                   \
    } while (0)

static int run_self_tests(void) {
    LOG_INFO("Running UK Vector Router self-tests...");

    /* 1. Transverse Mercator London -> Edinburgh distance check (~534 km chord) */
    float x1, y1, x2, y2;
    uk_project_tm(51.5074, -0.1278, &x1, &y1);
    uk_project_tm(55.9533, -3.1883, &x2, &y2);
    double tm_dist_km = sqrt((double)(x2 - x1) * (x2 - x1) + (double)(y2 - y1) * (y2 - y1)) / 1000.0;
    double hv_dist_km = haversine_m(51.5074, -0.1278, 55.9533, -3.1883) / 1000.0;
    TEST_CHK(fabs(tm_dist_km - hv_dist_km) < 2.0, "UK Transverse Mercator chord accuracy");

    /* 2. Seed UK demo in-memory SQLite, preprocess, verify contraction & peeling */
    sqlite3 *db = NULL;
    TEST_CHK(sqlite3_open(":memory:", &db) == SQLITE_OK, "sqlite :memory:");
    TEST_CHK(db_init_schema(db) == 0, "db_init_schema");
    TEST_CHK(db_seed_uk_demo(db) == 0, "db_seed_uk_demo");

    VmrRawGraph rg;
    raw_graph_init(&rg);
    TEST_CHK(db_load_raw_graph(db, &rg) == 0, "db_load_raw_graph");
    TEST_CHK(rg.node_count > 180, "raw UK demo nodes > 180");

    VmrGraph g;
    vmr_graph_init(&g);
    TEST_CHK(vmr_preprocess_graph(&rg, &g) == 0, "vmr_preprocess_graph");
    TEST_CHK(g.contracted_nodes > 120, "degree-2 chain contraction eliminated shape nodes");
    TEST_CHK(g.peeled_dead_nodes >= 10, "cul-de-sac 1-core peeling flagged dead ends");
    raw_graph_free(&rg);

    /* 3. Route London (1) -> Inverness (22) via Bidir, Unidir, Dijkstra */
    uint32_t s_lon = vmr_snap_coord(&g, 51.5074, -0.1278, 1);
    uint32_t t_inv = vmr_snap_coord(&g, 57.4778, -4.2247, 1);
    TEST_CHK(s_lon != VMR_INVALID && t_inv != VMR_INVALID, "spatial grid snap");

    VmrQuery qb = { MODE_VEC_BIDIR,  vmr_preset_profile("balanced"), 1 };
    VmrQuery qu = { MODE_VEC_UNIDIR, vmr_preset_profile("balanced"), 1 };
    VmrQuery qd = { MODE_DIJKSTRA,   vmr_preset_profile("balanced"), 1 };
    VmrResult rb, ru, rd;
    vmr_result_init(&rb);
    vmr_result_init(&ru);
    vmr_result_init(&rd);
    TEST_CHK(vmr_find_path(&g, s_lon, t_inv, &qb, &rb) == 0 && rb.found, "bidir London->Inverness");
    TEST_CHK(vmr_find_path(&g, s_lon, t_inv, &qu, &ru) == 0 && ru.found, "unidir London->Inverness");
    TEST_CHK(vmr_find_path(&g, s_lon, t_inv, &qd, &rd) == 0 && rd.found, "dijkstra London->Inverness");
    TEST_CHK(rb.cum_vec_cos_sim > 0.99, "cumulative vector sum reaches goal accurately");
    TEST_CHK(rb.nodes_expanded <= rd.nodes_expanded, "bidir expands <= dijkstra");
    vmr_result_free(&rb);
    vmr_result_free(&ru);
    vmr_result_free(&rd);

    /* 4. Verify O(1) island disconnection rejection (London -> Belfast) */
    uint32_t t_bel = vmr_snap_coord(&g, 54.5973, -5.9301, 1);
    VmrResult r_isl;
    vmr_result_init(&r_isl);
    TEST_CHK(vmr_find_path(&g, s_lon, t_bel, &qb, &r_isl) == 0, "island query returns cleanly");
    TEST_CHK(r_isl.found == 0 && r_isl.nodes_expanded == 0, "O(1) island rejection with 0 expansions");
    vmr_result_free(&r_isl);

    /* 5. Regression: a two-way dead-end stub must not decrement its live
     *    junction twice during 1-core peeling (old fwd+rev double decrement). */
    {
        VmrRawGraph tg;
        VmrGraph tgraph;
        raw_graph_init(&tg);
        TEST_CHK(raw_graph_add_node(&tg, 1, 51.0000, -1.0000) == 0, "peel node u");
        TEST_CHK(raw_graph_add_node(&tg, 2, 51.0010, -1.0000) == 0, "peel node v");
        TEST_CHK(raw_graph_add_node(&tg, 3, 51.0020, -1.0000) == 0, "peel node a");
        TEST_CHK(raw_graph_add_node(&tg, 4, 51.0010, -0.9990) == 0, "peel node b");
        TEST_CHK(raw_graph_add_node(&tg, 5, 51.0020, -0.9990) == 0, "peel node c");
        const float f_spd = 0.3f, f_cls = 0.3f, f_smo = 0.8f;
        TEST_CHK(raw_graph_add_edge(&tg, 0, 1, 100.0f, 12.0f, f_spd, f_cls, f_smo, RC_RESIDENTIAL, 0) == 0, "u->v");
        TEST_CHK(raw_graph_add_edge(&tg, 1, 0, 100.0f, 12.0f, f_spd, f_cls, f_smo, RC_RESIDENTIAL, 0) == 0, "v->u");
        TEST_CHK(raw_graph_add_edge(&tg, 3, 1, 100.0f, 12.0f, f_spd, f_cls, f_smo, RC_RESIDENTIAL, 0) == 0, "b->v");
        TEST_CHK(raw_graph_add_edge(&tg, 1, 2, 100.0f, 12.0f, f_spd, f_cls, f_smo, RC_RESIDENTIAL, 0) == 0, "v->a");
        TEST_CHK(raw_graph_add_edge(&tg, 2, 4, 100.0f, 12.0f, f_spd, f_cls, f_smo, RC_RESIDENTIAL, 0) == 0, "a->c");
        TEST_CHK(raw_graph_add_edge(&tg, 4, 2, 100.0f, 12.0f, f_spd, f_cls, f_smo, RC_RESIDENTIAL, 0) == 0, "c->a");
        TEST_CHK(raw_graph_add_edge(&tg, 3, 4, 100.0f, 12.0f, f_spd, f_cls, f_smo, RC_RESIDENTIAL, 0) == 0, "b->c");
        TEST_CHK(raw_graph_add_edge(&tg, 4, 3, 100.0f, 12.0f, f_spd, f_cls, f_smo, RC_RESIDENTIAL, 0) == 0, "c->b");

        vmr_graph_init(&tgraph);
        TEST_CHK(vmr_preprocess_graph(&tg, &tgraph) == 0, "peel regression preprocess");
        uint32_t ui = idmap_get(&tgraph.osm_map, 1);
        uint32_t vi = idmap_get(&tgraph.osm_map, 2);
        TEST_CHK(ui != VMR_INVALID && vi != VMR_INVALID, "peel regression node lookup");
        TEST_CHK((tgraph.node_flags[ui] & NODE_DEAD_TREE) != 0, "two-way stub peeled");
        TEST_CHK((tgraph.node_flags[vi] & NODE_DEAD_TREE) == 0, "live junction survives peeling");
        vmr_graph_free(&tgraph);
    }

    vmr_graph_free(&g);
    sqlite3_close(db);
    vmr_ws_free();
    LOG_INFO("All self-tests PASSED.");
    return 0;
}

static void print_usage(const char *prog) {
    printf(
        "UK Vector Router v%s (Preprocessed 8D SoA Vector-Corridor Pathfinding)\n"
        "Usage: %s [options]\n\n"
        "Data & Preprocessing:\n"
        "  --db <path>                   SQLite3 database path (default: uk_roadmap.db)\n"
        "  --init-uk-demo                Seed built-in UK inter-city demo map\n"
        "  --import-osm <file>           Import Geofabrik UK OSM XML (.osm) highway network\n"
        "  --import-pbf <file>           Import Geofabrik UK OSM PBF (.osm.pbf) highway network\n"
        "  --apply-turn-restrictions <file>  Add PBF turn restrictions to an existing DB\n"
        "  --import-csv <nodes> <edges>  Import CSV (id,lat,lon / u,v,class,kmh,smooth,flags)\n"
        "  --preprocess                  Run TM projection, contraction, 1-core peel, SoA pack\n"
        "  --save-bin <file.vgraph>      Save preprocessed SoA binary cache\n"
        "  --load-bin <file.vgraph>      Load preprocessed SoA binary cache (instant)\n\n"
        "Routing & Evaluation:\n"
        "  --route-coord <la1> <lo1> <la2> <lo2>   Snap WGS84 coords & route\n"
        "  --route-id <osmA> <osmB>                Route between OSM intersection IDs\n"
        "  --mode <bidir|unidir|astar|dijkstra|cch>\n"
        "                                        Search algorithm (default: bidir)\n"
        "  --profile <balanced|fastest|direct|highway|local>\n"
        "  --weights <a,s,c,t,m>                   Custom 5D weights (align,spd,cls,str,sm)\n"
        "  --benchmark                             Compare bidir vs unidir vs dijkstra\n"
        "  --geojson <out.geojson>                 Export computed route as GeoJSON\n"
        "  --interactive                           Start interactive UK routing shell\n"
        "  --self-test                             Run unit & integration tests\n",
        UVR_VERSION, prog);
}

int main(int argc, char **argv) {
    if (argc == 1) {
        print_usage(*argv);
        return 0;
    }

    const char *db_path    = "uk_roadmap.db";
    const char *osm_path   = NULL;
    const char *pbf_path   = NULL;
    int pbf_restrictions_only = 0;
    const char *csv_nodes  = NULL;
    const char *csv_edges  = NULL;
    const char *save_bin   = NULL;
    const char *load_bin   = NULL;
    const char *geojson_out= NULL;

    int do_demo = 0, do_preproc = 0, do_bench = 0, do_inter = 0, do_test = 0;
    int do_rcoord = 0, do_rid = 0;
    double la1 = 0, lo1 = 0, la2 = 0, lo2 = 0;
    int64_t id_a = 0, id_b = 0;
    SearchMode mode = MODE_VEC_BIDIR;
    VmrPreference pref = vmr_preset_profile("balanced");

    for (int i = 1; i < argc; i++) {
        const char *a = *(argv + i);
        if (strcmp(a, "--help") == 0 || strcmp(a, "-h") == 0) {
            print_usage(*argv); return 0;
        } else if (strcmp(a, "--self-test") == 0) {
            do_test = 1;
        } else if (strcmp(a, "--init-uk-demo") == 0) {
            do_demo = 1;
        } else if (strcmp(a, "--preprocess") == 0) {
            do_preproc = 1;
        } else if (strcmp(a, "--benchmark") == 0) {
            do_bench = 1;
        } else if (strcmp(a, "--interactive") == 0) {
            do_inter = 1;
        } else if (strcmp(a, "--db") == 0 && i + 1 < argc) {
            db_path = *(argv + (++i));
        } else if (strcmp(a, "--import-osm") == 0 && i + 1 < argc) {
            const char *p = *(argv + (++i));
            size_t pl = strlen(p);
            if (pl >= 4 && (strcmp(p + pl - 4, ".pbf") == 0 || strcmp(p + pl - 4, ".PBF") == 0)) {
                pbf_path = p;
            } else {
                osm_path = p;
            }
        } else if (strcmp(a, "--import-pbf") == 0 && i + 1 < argc) {
            pbf_path = *(argv + (++i));
        } else if (strcmp(a, "--apply-turn-restrictions") == 0 && i + 1 < argc) {
            pbf_path = *(argv + (++i));
            pbf_restrictions_only = 1;
        } else if (strcmp(a, "--import-csv") == 0 && i + 2 < argc) {
            csv_nodes = *(argv + (++i));
            csv_edges = *(argv + (++i));
        } else if (strcmp(a, "--save-bin") == 0 && i + 1 < argc) {
            save_bin = *(argv + (++i));
        } else if (strcmp(a, "--load-bin") == 0 && i + 1 < argc) {
            load_bin = *(argv + (++i));
        } else if (strcmp(a, "--geojson") == 0 && i + 1 < argc) {
            geojson_out = *(argv + (++i));
        } else if (strcmp(a, "--profile") == 0 && i + 1 < argc) {
            pref = vmr_preset_profile(*(argv + (++i)));
        } else if (strcmp(a, "--weights") == 0 && i + 1 < argc) {
            float wa, ws, wc, wt, wm;
            if (sscanf(*(argv + (++i)), "%f,%f,%f,%f,%f", &wa, &ws, &wc, &wt, &wm) == 5) {
                pref.w_align = wa; pref.w_speed = ws; pref.w_class = wc;
                pref.w_straight = wt; pref.w_smooth = wm;
                pref.pref_id = -1;
            }
        } else if (strcmp(a, "--mode") == 0 && i + 1 < argc) {
            const char *m = *(argv + (++i));
        if (c_stricmp(m, "bidir") == 0) mode = MODE_VEC_BIDIR;
        else if (c_stricmp(m, "unidir") == 0) mode = MODE_VEC_UNIDIR;
        else if (c_stricmp(m, "astar") == 0) mode = MODE_ASTAR;
        else if (c_stricmp(m, "dijkstra") == 0) mode = MODE_DIJKSTRA;
        else if (c_stricmp(m, "cch") == 0) mode = MODE_CCH;
        } else if (strcmp(a, "--route-coord") == 0 && i + 4 < argc) {
            do_rcoord = 1;
            la1 = strtod(*(argv + (++i)), NULL);
            lo1 = strtod(*(argv + (++i)), NULL);
            la2 = strtod(*(argv + (++i)), NULL);
            lo2 = strtod(*(argv + (++i)), NULL);
        } else if (strcmp(a, "--route-id") == 0 && i + 2 < argc) {
            do_rid = 1;
            id_a = strtoll(*(argv + (++i)), NULL, 10);
            id_b = strtoll(*(argv + (++i)), NULL, 10);
        } else {
            LOG_ERR("Unknown or incomplete option: %s", a);
            return 1;
        }
    }

    if (do_test && run_self_tests() != 0) {
        return 1;
    }

    int need_graph = (do_preproc || load_bin != NULL || do_rcoord || do_rid || do_inter || save_bin != NULL);
    int need_db    = (do_demo || osm_path != NULL || pbf_path != NULL || csv_nodes != NULL || (need_graph && load_bin == NULL));

    VmrGraph graph;
    vmr_graph_init(&graph);

    if (load_bin != NULL) {
        if (vmr_load_binary(&graph, load_bin) != 0) {
            return 1;
        }
    } else if (need_db) {
        sqlite3 *db = NULL;
        if (sqlite3_open(db_path, &db) != SQLITE_OK) {
            LOG_ERR("Cannot open SQLite '%s'", db_path);
            return 1;
        }
        if (db_init_schema(db) != 0) {
            sqlite3_close(db); return 1;
        }
        if (do_demo && db_seed_uk_demo(db) != 0) {
            sqlite3_close(db); return 1;
        }
        if (pbf_path && db_import_osm_pbf(db, pbf_path, pbf_restrictions_only) != 0) {
            sqlite3_close(db); return 1;
        }
        if (osm_path && db_import_osm_xml(db, osm_path) != 0) {
            sqlite3_close(db); return 1;
        }
        if (csv_nodes && csv_edges && db_import_csv(db, csv_nodes, csv_edges) != 0) {
            sqlite3_close(db); return 1;
        }
        if (do_preproc || do_rcoord || do_rid || do_inter || save_bin != NULL) {
            VmrRawGraph rg;
            raw_graph_init(&rg);
            if (db_load_raw_graph(db, &rg) != 0 || vmr_preprocess_graph(&rg, &graph) != 0) {
                raw_graph_free(&rg);
                sqlite3_close(db);
                vmr_graph_free(&graph);
                return 1;
            }
            raw_graph_free(&rg);
        }
        sqlite3_close(db);
    }

    if (save_bin != NULL && graph.node_count > 0) {
        vmr_save_binary(&graph, save_bin);
    }

    if (do_rcoord || do_rid) {
        uint32_t s = do_rcoord ? vmr_snap_coord(&graph, la1, lo1, 1)
                               : idmap_get(&graph.osm_map, id_a);
        uint32_t t = do_rcoord ? vmr_snap_coord(&graph, la2, lo2, 1)
                               : idmap_get(&graph.osm_map, id_b);
        if (s == VMR_INVALID || t == VMR_INVALID) {
            LOG_ERR("Failed to resolve start or goal location on UK graph");
            vmr_graph_free(&graph);
            return 1;
        }
        VmrQuery q = { mode, pref, 1 };
        VmrResult res;
        vmr_result_init(&res);
        if (vmr_find_path(&graph, s, t, &q, &res) == 0) {
            vmr_print_result(&graph, &res);
            if (geojson_out != NULL) {
                vmr_export_geojson(&graph, &res, geojson_out);
            }
            vmr_result_free(&res);
        }
        if (do_bench) {
            vmr_benchmark(&graph, s, t, &pref);
        }
    }

    if (do_inter && graph.node_count > 0) {
        run_interactive_shell(&graph);
    }

    vmr_ws_free();
    vmr_graph_free(&graph);
    if (g_mem.active != 0) {
        LOG_WARN("Unfreed allocations: %zu (%zu bytes)", g_mem.active, g_mem.cur_bytes);
    } else {
        LOG_DBG("Clean memory exit (peak %.2f MB)", (double)g_mem.peak_bytes / (1024.0 * 1024.0));
    }
    return 0;
}