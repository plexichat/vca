# VectorStar — Vector-Corridor A* (VCA*)

A single-file, dependency-free **VCA\*** router for large road networks.

VCA\* optimises a **direction-aware, preference-weighted vector cost** instead of
a single scalar such as travel time. Each edge carries a *direction* and a set of
*feature attributes* (road class, speed, surface). The cost of a path combines
how well each edge aligns with the remaining route direction (via two-hop
look-ahead continuation vectors) with a per-profile weighting of those
attributes. The result is a route that is straighter, keeps to preferred road
classes, and follows the overall bearing of the trip.

`algo.c` is one translation unit (C99, with optional OpenMP and AVX2 paths) that
does the whole pipeline — OSM XML / PBF / CSV / SQLite import, graph
contraction and preprocessing, and point-to-point search — with no dependencies
beyond `sqlite3.c` and `zlib`.

## Build

```sh
gcc -O2 -mavx2 -fopenmp -o uk_vector_router.exe algo.c sqlite3.o libz.a -I . -lm
# single-core (no OpenMP):
gcc -O2 -mavx2 -o uk_vector_router_sc.exe algo.c sqlite3.o libz.a -I . -lm
```

## Usage

```sh
# import + preprocess (PBF, OSM XML, or CSV) into a reusable .vgraph cache
./uk_vector_router.exe --import-pbf great-britain-latest.osm.pbf --db uk.db \
    --preprocess --save-bin uk.vgraph

# route between two coordinates
./uk_vector_router.exe --load-bin uk.vgraph \
    --route-coord 51.5074 -0.1278 55.9533 -3.1883 \
    --mode bidir --profile balanced

# eco route with custom fuel consumption (5.0 L/100km default)
./uk_vector_router.exe --load-bin uk.vgraph \
    --route-coord 51.5074 -0.1278 55.9533 -3.1883 \
    --profile eco --fuel 4.5

# scenic route (prefers winding, higher-class roads)
./uk_vector_router.exe --load-bin uk.vgraph \
    --route-coord 51.5074 -0.1278 55.9533 -3.1883 \
    --profile scenic

# shortest route (minimises distance, de-emphasises smoothness)
./uk_vector_router.exe --load-bin uk.vgraph \
    --route-coord 51.5074 -0.1278 55.9533 -3.1883 \
    --profile shortest

# route with live traffic multipliers (CSV: osm_id,multiplier)
./uk_vector_router.exe --load-bin uk.vgraph \
    --route-coord 51.5074 -0.1278 55.9533 -3.1883 \
    --traffic traffic.csv
```

**Modes** (`--mode`): `bidir` (default), `unidir`, `astar`, `dijkstra`, `cch`.

**Profiles** (`--profile`): `balanced` (default), `fastest`, `direct`, `highway`,
`local` (no-toll), `eco` (fuel-economy), `scenic` (winding roads), `shortest`
(minimise distance). Each profile sets the feature weight vector used in the cost.

**Fuel estimation:** `--fuel <l/100km>` sets the vehicle consumption rate (default
5.0 L/100km). The router estimates total fuel consumption based on distance and
road smoothness (lower smoothness = more stop-and-go = higher consumption).

**Live traffic:** `--traffic <file>` loads traffic multipliers from a CSV file
(format: `osm_id,multiplier` per line, `#` for comments). Multipliers are applied
per-edge to travel time, allowing the router to avoid congested roads.

**Environment overrides:** `VMR_ALT_SCALE` forces the ALT landmark heuristic
scale (default is chosen per query from the graph's admissible minimum); setting
`VMR_CCH=1` builds the CCH hierarchy at preprocess time.

## How the search works

- **Cost model.** Each edge stores a direction `(edge_ux, edge_uy)` and a scale
  `edge_attr_scale` derived from its feature attributes. The path cost weights
  edge length by a per-profile attribute score and rewards alignment between an
  edge and the two-hop continuation vector from the node it leads to.
- **Heuristic (ALT).** A small set of landmarks (31 on the UK graph) is chosen
  and their distances precomputed. The heuristic is the max over landmarks of
  `|d(L,t) - d(L,v)|` scaled by `attr_min`, which is admissible because every
  edge cost is at least `attr_min × length`.
- **Two-hop continuation.** For each node, a short forward/backward look-ahead
  vector is precomputed at preprocess time so the search can evaluate direction
  change without walking the path.
- **Contraction + reorder.** Low-degree nodes are contracted, dead-end components
  are peeled, and the graph is reordered by Hilbert curve so that each query
  touches a compact, cache-friendly region.

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

### All profiles (single core, London→Edinburgh)

| Profile | Dist | CPU | Notes |
|---|---|---|---|
| balanced | 597.75 km | 3–4 ms | Default; balances all features |
| fastest | 597.92 km | 3–4 ms | Prioritises speed |
| direct | 597.75 km | 3–4 ms | Prioritises straightness |
| highway | 597.92 km | 3–4 ms | Prioritises road class |
| local/no-toll | 598.49 km | 3–4 ms | Avoids tolls |
| **eco** | ~598 km | 3–4 ms | Fuel-economy; rewards smooth roads |
| **scenic** | ~600 km | 3–4 ms | Prefers winding, higher-class roads |
| **shortest** | ~590 km | 3–4 ms | Minimises distance |

**8 optimisation profiles, each < 20 ms per query (typically 3–4 ms).**

### Route quality (London→Edinburgh, balanced)

597.75 km (direct 533.5 km, **detour 1.12×**), 374 min, mean_cos **0.917**,
8-D alignment ~0.91, 919 hops.

## Comparison to incumbents

| | Dijkstra | A* (ALT) | CCH | **VCA\*** |
|---|---|---|---|---|
| Expanded (London→Edinburgh) | 1,239,374 | 926 | 2,877 | 1,644 |
| CPU (warm) | ~320 ms | ~1 ms | ~30 ms | **<1 ms** |
| Objective | time | time | time | **direction-aware vector** |
| Memory | — | — | large | **382 MB** |

VCA\* matches A\*'s speed on long routes and beats it on **route shape** (alignment,
straightness, class) at a fraction of CCH's memory.

## Status / caveats

- `.vgraph` format is v8; regenerate after changing the graph layout.
- CCH is opt-in (`VMR_CCH=1`) and experimental — the legacy CH builder does not
  scale cleanly to large graphs.
- The auto ALT scale is chosen per query from the graph's admissible minimum; for
  some short or mixed queries this is not the fastest admissible choice, so the
  search may expand more nodes than strictly necessary.
- Fuel estimation is a simple model (distance × rate × smoothness factor); it
  does not account for elevation, vehicle type, or driving style.
- Traffic multipliers are applied per-edge to travel time; they do not affect
  the ALT heuristic (which uses free-flow time).
- See `RESEARCH.md` for related work and the benchmark plan.