# Benchmarks

All measurements: AMD Ryzen 5 5600X (6c/12t), 32 GB RAM, Windows, MinGW-W64 GCC
`-O2 -mavx2 -fopenmp`, single core per query. Graph: UK `uk_roadmap.vgraph`
(2,714,206 nodes / 6,027,213 edges after contraction).

## Our router vs OSRM (industry standard), 24 inter-UK routes

OSRM = `router.project-osrm.org` public demo server (`driving` profile, live OSM).
"Our VCA" = bidirectional Vector-Corridor A* (`balanced` profile). The "ref" columns
are a time-optimal reference route (ALT A*) computed while benchmarking, shown only
to indicate how close VCA's vector-optimal route is to the shortest-time path.

| Route | OSRM km | ref km | VCA km | ref min | VCA min | Δdist | Δtime |
|---|---|---|---|---|---|---|---|
| London→Manchester | 321.0 | 319.6 | 312.7 | 191.3 | 201.7 | -2.6% | -17.3% |
| London→Edinburgh | 648.9 | 598.1 | 609.7 | 373.1 | 377.2 | -6.0% | -20.1% |
| London→Birmingham | 189.5 | 188.0 | 195.6 | 114.4 | 118.6 | +3.2% | -19.7% |
| Manchester→Edinburgh | 342.2 | 344.0 | 344.1 | 205.2 | 208.1 | +0.6% | -20.4% |
| Birmingham→Edinburgh | 460.8 | 462.2 | 463.0 | 266.3 | 270.7 | +0.5% | -19.8% |
| Bristol→Edinburgh | 591.8 | 593.4 | 593.4 | 337.2 | 337.2 | +0.3% | -20.8% |
| London→Cardiff | 241.5 | 240.9 | 244.1 | 147.3 | 163.0 | +1.1% | -12.9% |
| London→Newcastle | 455.2 | 443.8 | 454.8 | 258.7 | 277.1 | -0.1% | -15.5% |
| London→Liverpool | 340.2 | 338.8 | 333.0 | 198.7 | 211.1 | -2.1% | -16.6% |
| London→Leeds | 316.1 | 312.0 | 313.9 | 183.8 | 193.8 | -0.7% | -17.2% |
| Manchester→Birmingham | 139.5 | 138.6 | 138.4 | 82.1 | 88.9 | -0.7% | -16.0% |
| Manchester→Glasgow | 348.1 | 347.7 | 354.7 | 190.3 | 215.9 | +1.9% | -9.9% |
| Leeds→Edinburgh | 351.8 | 310.9 | 309.4 | 207.2 | 208.4 | -12.0% | -21.1% |
| Sheffield→London | 260.4 | 258.9 | 258.9 | 157.6 | 157.6 | -0.6% | -22.5% |
| Nottingham→London | 204.9 | 204.4 | 204.4 | 125.2 | 125.2 | -0.2% | -22.3% |
| Cambridge→Manchester | 290.0 | 258.6 | 262.9 | 161.0 | 163.6 | -9.4% | -21.8% |
| Oxford→Edinburgh | 578.8 | 580.0 | 580.0 | 331.2 | 331.2 | +0.2% | -21.0% |
| Norwich→London | 189.9 | 183.1 | 183.1 | 119.6 | 119.6 | -3.6% | -22.6% |
| Hull→London | 307.2 | 304.1 | 304.0 | 188.6 | 188.6 | -1.0% | -22.9% |
| Aberdeen→Edinburgh | 205.0 | 194.4 | 194.4 | 123.6 | 123.6 | -5.2% | -22.0% |
| Plymouth→London | 346.9 | 346.8 | 346.8 | 208.7 | 208.7 | -0.0% | -21.2% |
| Southampton→Edinburgh | 682.2 | 682.4 | 682.4 | 387.0 | 387.0 | +0.0% | -21.3% |
| Cardiff→Edinburgh | 623.0 | 624.4 | 624.4 | 359.4 | 359.4 | +0.2% | -21.0% |
| Liverpool→Glasgow | 354.8 | 354.4 | 352.3 | 196.4 | 196.7 | -0.7% | -20.5% |

**Mean VCA distance vs OSRM: −1.5%.** VCA routes are within ~1.5% of OSRM on
distance. The ~−20% time gap is because our `Est Time` is free-flow at the posted
speed, while OSRM models real-world turn/traffic costs. A few routes (Leeds→Edinburgh
−12%, Cambridge→Manchester −9.4%) diverge more because our graph is a specific UK
extract with a slightly different road set than OSRM's live planet data.

## Fuel-economy route (London→Edinburgh, car)

| Profile | Dist km | Time min | Fuel L | L/100km |
|---|---|---|---|---|
| balanced | 609.7 | 377.2 | 32.05 | 5.26 |
| **eco** | 594.3 | 406.2 | **31.44** | 5.29 |

The eco profile minimises a per-edge fuel model (U-shaped vs speed, smoothness
penalty) and genuinely finds a lower-fuel route (−1.9% vs balanced). It is slower to
compute (~0.32 s vs ~0.03 s) because the fuel-optimal path uses unavoidable slow
city edges; a fuel-Euclidean heuristic (`h = euclidean_km × min_fuel_per_km`) keeps
it admissible and cuts the expansions ~3.4× versus a time-only heuristic. Vehicle
type (speed cap + fuel curve) shifts both the route and the fuel figure — e.g.
`hgv-artic` (cap 85 km/h, 32 L/100km) yields a different path and ~210 L for the
same trip.

## Regeneration & live-update latency (single core unless noted)

The key point for live traffic: **traffic updates need no graph regeneration.**
Traffic is a query-time overlay (`traffic_mult`), so an update is just "reload the
file + re-route".

| Operation | Time | Notes |
|---|---|---|
| Load `.vgraph` cache (2.3 GB, mmap) | **194 ms** | cold; warm ~150 ms |
| Traffic file reload (10k lines) | **< 5 ms** | O(lines) via osm_id index; 24 MB array reused |
| Re-route, vca-bidir | **~27 ms** | warm |
| **Live traffic update cycle** (reload + re-route) | **~30–50 ms** | no regeneration |
| Full graph regen: SQLite DB load (1.8 GB, 17 M nodes) | ~38 s | I/O-bound |
| Full graph regen: preprocess (OpenMP 6 threads) | ~10 s | contraction + reorder + two-hop + ALT |
| Full graph regen: save `.vgraph` (1.66 GB) | ~1–2 s | |
| **Full graph regen total** (DB → .vgraph) | **~50 s** | ~110 s single-core |

So a live-traffic update lands in **~30–50 ms** (file reload + re-route), versus
**~50 s** if the road network itself changes and the graph must be rebuilt. The
one-time graph build is dominated by reading the 17 M-node SQLite DB; the actual
routing preprocess is only ~10 s with OpenMP.