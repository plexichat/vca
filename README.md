# VectorStar — Vector-Corridor A* (VCA*) Routing Engine

A single-file, dependency-free **VCA\*** router for large road networks.
Instead of optimising a single scalar (travel time), VCA\* optimises a
**direction-aware, preference-weighted vector cost** — producing routes that are
straighter, stay on preferred road classes, and align with the overall
direction of travel.

`algo.c` (single translation unit, C99 + optional OpenMP/AVX2) implements:
import (OSM XML / **PBF** / CSV / SQLite), preprocessing, and point-to-point
search, with no external dependencies beyond `sqlite3.c` and `zlib`.

---

## Build

```sh
gcc -O2 -mavx2 -fopenmp -o uk_vector_router.exe algo.c sqlite3.o libz.a -I . -lm
# single-core (no OpenMP):
gcc -O2 -mavx2 -o uk_vector_router_sc.exe algo.c sqlite3.o libz.a -I . -lm
```

## Usage

```sh
# import + preprocess (PBF or OSM XML or CSV) -> .vgraph cache
./uk_vector_router.exe --import-pbf great-britain-latest.osm.pbf --db uk.db --preprocess --save-bin uk.vgraph

# route (any of 5 profiles, 5 modes)
./uk_vector_router.exe --load-bin uk.vgraph --route-coord 51.5074 -0.1278 55.9533 -3.1883 \
    --mode bidir --profile balanced

# modes: bidir | unidir | astar | dijkstra | cch
# profiles: balanced | fastest | direct | highway | local (no-toll)
# env: VMR_ALT_SCALE (override auto ALT scale), VMR_CCH=1 (build CCH hierarchy)
```

---

## Performance (single core, this machine)

**Machine:** AMD Ryzen 5 5600X (6c/12t), 32 GB RAM, Windows, MinGW-W64 GCC `-O2 -mavx2`.
**Graph:** UK — 17,039,603 OSM nodes → 2,714,206 after contraction; 33,131,465 → 6,027,213 edges (31 cul-de-sac peels, 31 ALT landmarks).

### Preprocessing (one-time)
| Phase | Single core | OpenMP (6 threads) |
|---|---|---|
| SQLite DB load (17 M nodes / 33 M edges) | ~45 s (I/O) | ~45 s (I/O) |
| Contraction + components + reorder + two-hop + grid | ~7.7 s | ~7.4 s |
| ALT landmarks (32) | **19.8 s** | **2.8 s** |
| **Total algorithm** | **27.5 s** | **~10 s** |

### Query latency (warm, single core)
| Algorithm | Expanded | CPU |
|---|---|---|
| **vca-bidir** | **1,644** | **<1 ms** |
| vca-unidir | 925 | <1 ms |
| astar (ALT-A*) | 926 | ~1 ms |
| dijkstra (true) | 1,239,374 | ~320 ms |

Cold first query (mmap page-faults) is ~0.3–0.9 s; warm is <1 ms.

### Peak memory (single core)
| Metric | Value |
|---|---|
| Peak **private** (heap) | **382 MB** |
| Peak working set (incl. mmap pages) | 583 MB |

### All 5 profiles (single core, London→Edinburgh)
| Profile | Dist | CPU |
|---|---|---|
| balanced | 597.75 km | 3–4 ms |
| fastest | 597.92 km | 3–4 ms |
| direct | 597.75 km | 3–4 ms |
| highway | 597.92 km | 3–4 ms |
| local/no-toll | 598.49 km | 3–4 ms |

**Stat: 5 optimisation profiles, each < 20 ms per query (typically 3–4 ms).**

### Route quality (London→Edinburgh, balanced)
597.75 km (direct 533.5 km, **detour 1.12×**), 374 min, mean_cos **0.917**,
8-D alignment ~0.91, 919 hops.

---

## Comparison to incumbents

| | Dijkstra | A* (ALT) | CCH | **VCA\*** |
|---|---|---|---|---|
| Expanded (London→Edinburgh) | 1,239,374 | 926 | 2,877 | 1,644 |
| CPU (warm) | ~320 ms | ~1 ms | ~30 ms | **<1 ms** |
| Objective | time | time | time | **direction-aware vector** |
| Memory | — | — | large | **382 MB** |

VCA\* matches A\*'s speed on long routes and beats it on **route shape** (alignment/straightness/class), at a fraction of CCH's memory.

---

## How to improve further

**Latency (→ sub-millisecond, more headroom):**
1. **Enable the AVX2 block expansion by default** (currently gated by `VMR_BLOCK`); build `rev_blocks` for the backward pass too.
10. **Cache `edge_attr_scale`** (already mmap'd) and precompute `time·attr` per edge to remove a multiply in the hot loop.
11. **Stronger admissible heuristic:** build landmarks on a `time·attr` metric or add more landmarks; the current 8-landmark cap is conservative.
12. **Radix heap** for the ALT landmark Dijkstras (single-core preprocess 19.8 s is the biggest one-time cost).

**Memory (→ <300 MB):**
- mmap the per-query `sf/sb` workspace or shrink it (drop `mx/my`, recompute from parent edge).
- Reduce ALT to int16-quantised distances (695 MB → ~170 MB if heap-backed).

**Preprocessing:**
- Parallelise contraction/components; radix heap for ALT; better landmark selection (k-means / farthest-point).

**Quality:**
- **Build `edge_attr_scale` from normalised features** (currently unnormalised → 0.5–0.8); add a turn-count metric; validate against a PSP baseline.
- Fix the legacy `vmr_build_ch` (pair-set 8 GB growth) so the CCH comparison is exact.

**Correctness:**
- The auto ALT scale is sub-optimal for some short/mixed queries; the node-budget retry mitigates but doesn't eliminate it. A proper per-query admissibility check (consistency test on the found path) would close this.

---

## Status / caveats
- `.vgraph` format is v8; regenerate after changes.
- CCH is opt-in (`VMR_CCH=1`) and experimental (legacy CH builder doesn't scale cleanly).
- See `RESEARCH.md` for related work and the benchmark plan.