# Vector-Corridor Routing: A Direction-Aware, Preference-Weighted Shortest-Path Engine

*A research note / paper skeleton. Working title: "Vector-Corridor Routing: Preference-Weighted, Direction-Aware Shortest Paths at CH-Class Latency."*

---

## Abstract

We present a point-to-point router for large road networks that optimises a
**direction-aware, preference-weighted vector cost** rather than a single scalar
(travel time). Each edge carries an 8-D vector `(time, dx, dy, ux, uy, f_speed,
f_class, f_smooth)`; the cost combines a 5-D user preference (distance, speed,
road-class, straightness, smoothness) with **alignment** (cosine to the goal),
**directional momentum** (heading continuity), and a precomputed **2-hop
continuation cosine** (does the edge continue toward the goal over the next
edge). The search is a **bidirectional vector-corridor A\*** whose ALT landmark
heuristic is **auto-scaled per query** to hold latency under 20 ms. On the UK
(c. 17 M OSM nodes) it answers London→Edinburgh in <1 ms at ~384 MB peak RSS;
all five preference profiles answer in 3–4 ms. We position the contribution
against prior art and argue the novelty is an *engineering synthesis* (PSP +
turn/straightness costs + ALT) with two genuinely new ingredients: the 2-hop
continuation cosine as a **cost term** and the admissibility-aware ALT
auto-scaling for a non-metric objective.

---

## 1. Motivation

Commercial routers optimise a vendor-defined scalar (usually travel time). This
is fast but produces routes that ignore *shape*: they zig-zag, change road class
frequently, and do not correspond to how people describe or prefer routes. Prior
work on "natural"/"fewest-turn"/"simplest" routing addresses this but at the cost
of much slower queries or metric-specific preprocessing.

Goal: produce **straighter, class-aware, direction-consistent** routes at
**CH-class latency (<20 ms)** and **low memory (~400 MB)**, with the ability to
change the preference at query time (no rebuild).

---

## 2. Related work (prior art)

| Mechanism in this engine | Established prior art | Reference |
|---|---|---|
| Per-edge cost **vector** + user preference **linear combination** | **Personalized Shortest Path (PSP)** / parametrized costs | Dibbelt, Strasser, Wagner, *Fast Exact Shortest Path … with Parametrized Costs*, SEA 2015 (arXiv:1509.03165) |
| Multi-criteria Pareto / regret filtering | Multi-criteria route planning | *Multi-Criteria Route Planning with Little Regret*, ATMOS 2025 |
| Arbitrary metric + turn costs + fast queries | **Customizable Route Planning (CRP)** / CCH | Delling, Goldberg, Pajor, Werneck, 2011–2013 (Bing Maps) |
| Straightness / fewest-turn / "natural roads" cost | Fewest-turn / simplest-path routing | Jiang & Liu, *Fewest-turn map directions based on natural roads*, IJGIS 2011 (arXiv:1003.3536); *Routing Directions: Keeping it Fast and Simple*, 2021 |
| ALT landmark heuristic + bidirectional | ALT / CH / CRP speedups | Goldberg & Harrelson, ALT 2005; *Optimizing Landmark-Based Routing*, 2013 |
| 2-hop look-ahead | k-step look-ahead A* heuristic; 2-hop labeling | Chen, *Improved A\* … k-step Look-ahead*, 2021; H2H; QHL, SIGMOD 2023 |
| Learned / custom cost functions | Learned routing costs | *Customizable and Jointly Optimized Route Planning*, 2025/26 |
| "Vector field toward goal" | Navigation fields / vector-field planners (robotics) | Mesh Navigation (CVP); Patil et al., *Navigation Fields*, 2010; Boroujeni et al., 2018 |

**Takeaway:** every *ingredient* is known. The cost is a **Personalized Shortest
Path** instance; the straightness/continuity terms are **fewest-turn / natural
routes**; the speedups are **ALT + bidirectional**; turn handling is **CRP-style**.

---

## 3. What is (and isn't) new

### Not new
- Preference-weighted linear cost → **PSP**.
- Straightness / fewest-turns / natural-road continuity → **Jiang & Liu 2010**.
- ALT + bidirectional A*, turn costs, customizable metrics → **CRP/CH/CCH**.

### Genuinely new (this engine)
1. **2-hop continuation cosine as a *cost* term.** Precomputed per edge
   (`th_sx/th_sy/th_inv`): the best alignment over all 2-edge continuations.
   Prior look-ahead work uses look-ahead to *tighten a heuristic*; here it is a
   *routing objective* (smoother, more "natural" paths).
2. **Admissibility-aware ALT auto-scaling for a non-metric cost.** The 8-D cost
   is not a metric, so a time-based ALT must be scaled by `min(attr_scale)` to
   stay admissible. We auto-select the scale per query (route-length threshold +
   node-budget retry) to keep latency <20 ms without sacrificing correctness.
7. **A single-file, dependency-free SoA + AVX2 + mmap implementation** at
   ~384 MB / <20 ms (an *artifact/systems* contribution).

### Honest verdict
This is an **incremental engineering synthesis**, not a new algorithmic paradigm.
The most defensible *research* claim is **route quality**: alignment/straightness
vs. PSP at CH-class latency. To be publishable, the 2-hop-continuation cost and
the non-metric ALT auto-scaling must be **quantitatively compared** to PSP/CRP
baselines (see §5).

---

## 4. Approach (summary)

1. **8-D SoA graph**: per-edge `(time, dx, dy, ux, uy, f_speed, f_class,
   f_smooth)` + 5 profile `edge_attr_scale` arrays.
7. **Cost** = `time · attr_scale · align_pen(cos1,cos2,mom) · corridor_pen`,
   with the 5-D preference folded into `attr_scale`.
8. **2-hop continuation cosine** precomputed into `th_*` / `rth_*` tables and
   persisted in the `.vgraph` (mmap'd at load, no rebuild).
7. **Bidirectional vector-corridor A\*** with an **ALT** heuristic; the ALT
   scale is auto-selected (strong for long/motorway routes, admissible
   `min(attr)` for short/mixed, node-budget fallback).
8. **mmap** for ALT / two-hop / CH; compact 8-landmark AVX2 ALT (64 B/relax).

---

## 5. Benchmark plan (to make it publishable)

**Baselines** (not just Dijkstra):
- Dijkstra, A* (Euclidean), **ALT** (time), **CRP/CCH** (time), **PSP** (5-D
  linear cost, via Dijkstra on α·c), and the "fewest-turn / natural roads"
  objective.

**Datasets**: UK (done), US (in progress), plus a few European regions.

**Metrics**:
- Latency (warm/cold), peak RSS.
- **Route quality**: detour ratio, mean alignment cosine, 8-D alignment,
  **turn count**, road-class mix, vs. the time-optimal path.
- Admissibility: fraction of queries where the heuristic over-estimates
   (sub-optimality rate) as a function of the auto-scale.

**Ablations**: 2-hop on/off; ALT scale (0.5→1.0); landmarks 8/16/32; mmap vs heap.

**Headline claim to prove**: *direction-aware 2-hop-continuation cost yields
straighter/fewer-turn routes than PSP at comparable or better latency.*

---

## 6. Key stats (measured, UK, Ryzen 5 5600X, 32 GB)

- **Profiles (optimization classes): 5** — `balanced, fastest, direct, highway,
  local/no-toll` — each answers London→Edinburgh in **3–4 ms** (<20 ms).
- **Road classes weighted: 6** — motorway, trunk, primary, secondary, tertiary,
  minor.
- **Latency**: vec-bidir London→Edinburgh **<1–5 ms** (1,644 expansions);
  Cardiff 19 ms; Manchester <1 ms (auto-scale).
- **Peak RSS**: **382 MB private / 657 MB working set** (mmap'd ALT+two-hop+CH).
- **Load**: `.vgraph` 1.3 GB (UK) mmap-load in 0.17 s; two-hop rebuild 0.4 s.
- **vs baselines**: true Dijkstra 1.24M expansions / ~250 ms; ALT-A\* 926 exp /
  ~1 ms; CCH 2,877 exp / ~30 ms.

---

## 7. Reproduce

```
# build
gcc -O2 -mavx2 -o uk_vector_router.exe algo.c sqlite3.o libz.a -I . -lm

# preprocess UK (or US)
./uk_vector_router.exe --db data/europe/uk/uk_roadmap.db --preprocess --save-bin uk.vgraph

# route (any of 5 profiles, 5 modes)
./uk_vector_router.exe --load-bin uk.vgraph --route-coord 51.5074 -0.1278 55.9533 -3.1883 --mode bidir --profile balanced --benchmark

# env knobs
VMR_ALT_SCALE=0.88   # override auto ALT scale
VMR_CCH=1            # build the CCH comparison hierarchy
```

---

## 8. Limitations / threats to validity

- The strong-scale auto-select is **sub-optimal for some short/mixed queries**
  (Cardiff); the node-budget retry mitigates but does not eliminate this.
- `vmr_build_ch` (legacy CH) does not scale cleanly to 16 M+ edges (8 GB
  reallocs during pair-set growth) — needs a rewrite for the CCH baseline.
- Quality claims need a turn-count metric and a PSP baseline (not yet measured).