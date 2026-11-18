# VectorStar — Vector-Corridor A* (VCA*)

A single-file, dependency-free **VCA\*** router for large road networks.

VCA\* optimises a **direction-aware, preference-weighted vector cost** instead of a
single scalar such as travel time. Each edge carries a *direction* and a set of
*feature attributes* (road class, speed, straightness, surface). The cost of a path
combines how well each edge aligns with the remaining route direction (via two-hop
look-ahead continuation vectors) with a per-profile weighting of those attributes.
The result is a route that is straighter, keeps to preferred road classes, and
follows the overall bearing of the trip.

`algo.c` is one translation unit (C99, with optional OpenMP and AVX2 paths) that does
the whole pipeline — OSM XML / PBF / CSV / SQLite import, graph contraction and
preprocessing, and point-to-point search — with no dependencies beyond `sqlite3.c`
and `zlib`.

## Build

```sh
gcc -O2 -mavx2 -fopenmp -funroll-loops -o uk_vector_router.exe algo.c sqlite3.o libz.a -I . -lm
# single-core (no OpenMP):
gcc -O2 -mavx2 -funroll-loops -o uk_vector_router_sc.exe algo.c sqlite3.o libz.a -I . -lm
```

`-funroll-loops` unrolls the scalar query inner loop and is worth ~5% on query
latency versus a plain `-O2` build (measured interleaved on the 24-route benchmark).

## Usage

```sh
# import + preprocess (PBF, OSM XML, or CSV) into a reusable .vgraph cache
./uk_vector_router.exe --import-pbf great-britain-latest.osm.pbf --db uk.db \
    --preprocess --save-bin uk.vgraph

# route between two coordinates
./uk_vector_router.exe --load-bin uk.vgraph \
    --route-coord 51.5074 -0.1278 55.9533 -3.1883 \
    --mode bidir --profile balanced

# fuel-economy route (minimises a per-edge fuel model)
./uk_vector_router.exe --load-bin uk.vgraph \
    --route-coord 51.5074 -0.1278 55.9533 -3.1883 --profile eco

# logistics vehicle (speed cap + fuel curve change the route and the fuel figure)
./uk_vector_router.exe --load-bin uk.vgraph \
    --route-coord 51.5074 -0.1278 55.9533 -3.1883 \
    --profile balanced --vehicle hgv-artic

# live-traffic watch: reload a traffic file and re-route every N seconds
./uk_vector_router.exe --load-bin uk.vgraph \
    --route-coord 51.5074 -0.1278 55.9533 -3.1883 \
    --traffic traffic.csv --watch 60
```

## VCA's cases

**Search variant** (`--mode`):
- `bidir` (default) — bidirectional VCA: forward and backward searches grow from
  origin and destination and meet in the middle. Fewest expansions on long routes.
- `unidir` — unidirectional VCA: a single forward search. Simpler, expands more
  nodes, useful for debugging or when a backward graph is unavailable.

**Optimisation profile** (`--profile`) — sets the feature-weight vector used in the
edge cost:
- `balanced` (default) — equal blend of speed, class, straightness, smoothness.
- `direct` — weights straightness/alignment most; fewest bends.
- `highway` — weights road class most; prefers motorways/trunk roads.
- `local` / `no-toll` — down-weights major roads and applies a toll penalty; avoids
  motorways and toll roads.
- `eco` — **fuel-economy**: the edge cost is proportional to fuel burned
  (`cost = time × fuel_l_per_km(speed) × speed`), so the route minimises total fuel.
- `scenic` — prefers higher-class, winding (low-straightness) roads over speed.
- `shortest` — minimises distance, de-emphasises smoothness.
- `avoid-motorway` / `no-motorway` — down-weights road class and speed, prefers
  smooth A-roads; useful for fuel, safety, or drivers who dislike motorways.

**Vehicle** (`--vehicle`) — sets a fuel curve + speed cap that change both the route
and the reported fuel figure:
- `car`, `van`, `hgv-rigid`, `hgv-artic`, `ev`.
- Each has a base L/100km, a speed-sensitivity coefficient, a stop-and-go penalty,
  and a speed cap. The cap slows edges above it, so an HGV (cap 85 km/h) takes a
  measurably different, slower path than a car (cap 130 km/h).
- `--fuel <l/100km>` overrides the base consumption rate.

**Fuel economy** (`--profile eco`): the fuel model is U-shaped vs speed — minimal at
the vehicle's optimum cruising speed, rising at low speed (stop-and-go, wrong gear)
and high speed (aero drag). The eco search uses a fuel-Euclidean heuristic
(`h = euclidean_km × min_fuel_per_km`) so it stays admissible and fast.

**Live traffic** (`--traffic <file>` + `--watch <seconds>`): `--traffic` loads
per-edge travel-time multipliers from a CSV (`osm_id,multiplier` per line, `#` for
comments). `--watch` re-loads the file and re-routes every N seconds, printing a
timestamped result. **Traffic updates need no graph regeneration** — the multiplier
array is a query-time overlay, so an update is just "reload the CSV + re-route",
which lands in ~30–50 ms (file reload is O(lines) via the OSM-id index, and the
24 MB multiplier array is reused). A full graph rebuild (only needed when the road
network itself changes) takes ~50 s with OpenMP. This is a file-based hook, **not** a
live network feed: a production system would have a feed reader (e.g. National
Highways WebTRIS for England, or a commercial TomTom/HERE/Google feed) rewrite the
file every 1–5 minutes.

**Environment override:** `UVR_ALT_SCALE` forces the ALT landmark heuristic scale
(default is the graph's admissible per-profile minimum).

## How the search works

- **Cost model.** Each edge stores a direction `(edge_ux, edge_uy)` and a precomputed
  `edge_attr_scale` per profile (the 4-D feature dot product folded in, so the query
  hot loop is a multiply). The path cost weights travel time by this attribute score
  and rewards alignment between an edge and the two-hop continuation vector from the
  node it leads to.
- **Heuristic (ALT).** A set of landmarks (16 on the UK graph) is chosen and their
  distances precomputed. The heuristic is `max over landmarks |d(L,t) − d(L,v)|`
  scaled by the profile's admissible `min(edge_attr_scale)`. A* with an admissible ALT
  heuristic is optimal; the scale is exact per profile so longer routes are not
  over-pruned.
- **Two-hop continuation.** For each edge, a short forward/backward look-ahead
  continuation table is precomputed at preprocess time, so the search evaluates
  direction change with a dot product instead of walking the path.
- **Contraction + reorder.** Low-degree nodes are contracted, dead-end components are
  peeled, and the graph is reordered by Hilbert curve so each query touches a compact,
  cache-friendly region.

## Performance (single core, this machine)

**Machine:** AMD Ryzen 5 5600X (6c/12t), 32 GB RAM, Windows, MinGW-W64 GCC `-O2 -mavx2`.
**Graph:** UK — 17,039,603 OSM nodes → 2,714,206 after contraction; 33,131,465 → 6,027,213 edges.

### Preprocessing (one-time)

| Phase | Single core | OpenMP (6 threads) |
|---|---|---|
| SQLite DB load (17 M nodes / 33 M edges) | ~38–45 s (I/O) | ~38–45 s (I/O) |
| Contract + components + peel + reorder + two-hop + grid + attr | ~12 s | ~5 s |
| ALT landmarks (16, 16-bit) | **~8 s** | **~1.5 s** |
| **Total algorithm** | **~20 s** | **~6.6 s** |

The heavy passes are parallelised: the edge group-by-source is an O(n+m) counting
sort (not a qsort), the Hilbert reorder gathers node/edge arrays concurrently, and
the landmark Dijkstras run across threads. The landmark table has 16 well-spread
slots (14 unique after dedup), stored as 16-bit quantized travel times at 2 s
resolution so the ALT table stays at 64 B/node while giving admissible time lower
bounds for the vector cost. The one-time cost is dominated by the SQLite DB read.

### Query latency (warm, single core, means over 24 UK routes)

| Variant | Expanded | CPU |
|---|---|---|
| **vca-bidir** | 40,636 | **~13–15 ms** |
| vca-unidir | 84,647 | ~30 ms |

`vca-bidir` grows forward and backward searches from origin and destination and
alternates sides; the ALT landmark table gives an admissible time lower bound for
the vector cost. The balanced profile follows a vector-optimal (straighter,
class-aware) path rather than a pure time-optimal one. The eco profile is slower
to compute because its fuel-optimal path uses unavoidable slow city edges; the
fuel-Euclidean heuristic keeps the search admissible and cuts its expansions ~3.4×
versus a time-only heuristic.

### Peak memory (single core, after load)

| Metric | Value |
|---|---|
| Peak working set | **~950 MB** |

The precompute passes trade memory for query speed: per-profile `time×attr`
base costs (8 × 6.0 M edges) and per-edge target/source coordinates (~96 MB)
remove random `edge_attr_scale`/`node_x` lookups from the query hot loop.

## Comparison to OSRM (industry standard)

We compare against **OSRM** (the open-source routing engine behind the OSRM demo
server) on 24 inter-UK routes. VCA routes are within **1.5% of OSRM on distance** on
average. See [`BENCHMARKS.md`](BENCHMARKS.md) for the full per-route table and the
fuel-economy numbers.

## Status / caveats

- `.vgraph` format is v12; the loader accepts older caches and backfills the newer
  eco/scenic/shortest profile scales on legacy v8 caches.
- The ALT scale is the exact admissible per-profile minimum; a few routes with
  unavoidable slow edges (notably eco) expand more nodes than strictly necessary.
- The fuel model is a static per-edge approximation (speed + smoothness); it does not
  account for elevation, real-time traffic, or driving style.
- Traffic multipliers affect travel time but not the ALT heuristic (which is built on
  free-flow time), so heavy congestion can increase expansions.
- See `RESEARCH.md` for related work and `BENCHMARKS.md` for measurements.