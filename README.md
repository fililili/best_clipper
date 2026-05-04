# best_clipper

A C++20 header-only library for polygon boolean operations (union, intersection, XOR, etc.) based on **snap rounding** and **winding number face traversal**.

## Core Ideas

### 1. Snap Rounding with Spatial Index

Snap rounding is a technique to convert arbitrary-precision geometric arrangements into integer-grid representations while preserving topological correctness. The theoretical foundation comes from:

> **Guibas, L. J., & Marimont, D. H. (1998).** *Rounding Arrangements Dynamically.* International Journal of Computational Geometry & Applications, Vol. 8, No. 2, pp. 157-176.

The key idea: all segment endpoints and intersection points are treated as **hot pixels** (integer grid points). Every segment is then "snapped" to pass through these hot pixels, producing a planar subdivision that is topologically equivalent to the true arrangement. A valid snap rounding guarantees:

- No two segments cross between grid lines without an intermediate hot pixel
- No segment passes closer than half a pixel to a hot pixel without passing through it
- The resulting planar graph is a valid subdivision of the plane

Unlike the original paper which handles dynamic (online) insertion, this library processes all segments in batch (offline). Instead of a dynamic rounding structure, a **spatial index tree** (Boost.Geometry rtree) is used to accelerate segment-intersection and point-on-segment queries. All geometric computations use **strict integer arithmetic** — no floating point inaccuracy.

Pipeline for snap rounding:
1. Find all segment-segment intersection points, snapped to integer grid (hot pixels)
2. For each original segment, find all hot pixels lying on it
3. Split each segment into sub-segments between consecutive hot pixels
4. Deduplicate edges (normalized so start < end)

### 2. Winding Number

**Definition.** For a point P in the plane, its **winding number** w(P) is defined as: shoot a ray from P to infinity in any direction (avoiding vertices), count each crossing of an input polygon's directed boundary — **+1** when the boundary crosses left-to-right, **-1** when right-to-left. The sum is w(P).

A face is a maximal connected region of constant winding number. Therefore:

| Operation | Winding Number Condition |
|-----------|-------------------------|
| OR (union) | w > 0 |
| AND (intersection) | w > 1 |
| XOR | w == 1 |
| A - B | w_A > 0 and w_B = 0 |

#### Edge Power

Each edge has a **power** computed from its direction relative to the hot pixel index order. When `start < end` (edge traverses from lower to higher pixel index), power = +1; when reversed, power = -1. Power represents the winding number jump when crossing the edge from right face to left face: `w(left) - w(right) = -power`. When edges from different input polygons overlap, their powers are **summed** (cancellation when equal and opposite). This is computed during snap rounding deduplication.

#### Chain Building

Adjacent edges with the **same power** and **no branching** (each intermediate vertex has exactly one incoming and one outgoing edge with the same power) are collapsed into a **chain**. This significantly reduces the graph size for subsequent steps.

#### Half Chains (Half-Edge Data Structure)

Each chain is split into two **half chains**, which form a standard **half-edge data structure**:

- **Forward half chain** (even index): traverses the chain in its original direction, belonging to the face on its **right** side
- **Reverse half chain** (odd index): traverses the chain in the opposite direction, belonging to the face on its **right** side (the other face of the chain)

This matches **Boost.Geometry's convention**: every half-edge belongs to the face on its right. (A valid Boost outer ring — math-CW — has the interior on the right of each edge.)

**Why the half chain convention must match the input/output convention.** The input polygon segments carry directional information: edge power is derived from the traversal direction, which for a Boost-valid ring already encodes which side is interior. If the half chain convention differed (e.g. left-face), then following `next_half_chain` around a face would trace rings in the **opposite** direction from the input, producing output rings with reversed orientation — necessitating a `std::reverse` on every ring. By aligning the half chain convention with Boost's convention:
- Face traversal via `next_half_chain` produces rings with the **same orientation** as the input
- Output rings are Boost-valid **without any post-processing**
- Edge power, `next_half_chain` traversal, and ring orientation all share a consistent directional framework

The key operations:

| Half-Edge Concept | Implementation |
|------------------|----------------|
| twin / opposite | `half_chain.dual()` = `half_chain.id ^ 1` |
| next (around face) | `next_half_chain[half_chain.id]` = next half-edge in the same face |
| power (winding delta) | forward: `+chain_power`, reverse: `-chain_power` |

#### Coplanarity (Same Face)

Two half chains are **coplanar** iff they belong to the same face (same winding number). Coplanarity comes from three sources:

1. **Sector coplanarity**: At each vertex, departing half chains are sorted by CCW angle. The sector between adjacent half chains `prev` and `cur` (CCW) lies on the right of `cur`. Since a half chain belongs to its right face, `cur` and `prev.dual()` belong to the same face: `coplanar(cur, prev.dual())`.

2. **Ray-casting coplanarity**: From the leftmost vertex of each component, cast a ray in the -x direction (with an infinitesimal +dy perturbation to avoid passing through vertices). The ray is a continuous curve. Since it does not cross any boundary until it hits something, every point along the ray belongs to the **same face**. Therefore the half chain whose right face contains the ray's origin and the half chain on the hit side (the side the ray approaches from) are coplanar. If the ray hits nothing, that half chain's right face is the exterior (winding = 0). This seeds winding number propagation.

3. **Dual cancellation coplanarity**: When both forward and reverse half chains of a chain survive filtering, they cancel each other. Their adjacent faces merge: `coplanar(forward, reverse)`.

#### Winding Number Propagation

A weighted graph is built on half-chain IDs with three edge types:

| Edge type | Weight | Meaning |
|-----------|--------|---------|
| coplanar pair | 0 | same winding number |
| ray coplanar pair | 0 | same winding number (exterior linkage) |
| dual (`i` ↔ `i^1`) | `-power(i)` | `winding[i^1] = winding[i] - power(i)` |

A **DFS** from the exterior seeds (identified by ray casting, `winding = 0`) propagates winding numbers to all reachable half chains. Each half chain's winding number is the winding of its right face.

#### Output Construction

1. Filter half chains: keep those whose `winding` satisfies the operation (e.g. `winding > 0` for union)
2. Dual cancellation: if both forward and reverse half chains of a chain survive, kill both, merge faces via `coplanar(forward, reverse)`
3. Rebuild `next_half_chain` respecting surviving half chains and face boundaries
4. Trace surviving half chains along `next_half_chain` into closed rings
5. Assemble rings into polygons (max-area ring = outer, rest = holes)

#### Self-Intersecting Input

Self-intersecting input polygons are handled naturally: the snap rounding step splits all intersecting segments at their intersection points, creating a proper planar subdivision. The winding number of each resulting face correctly reflects the net enclosure count. Use `self_or()` to first compute the union of a self-intersecting polygon with itself, which resolves self-intersections into a valid multi-polygon.

## Algorithm Pipeline Summary

```
Input polygons
  → Collect segments
  → Snap rounding (find intersections, split segments, deduplicate edges)
  → Assign edge powers (+1/-1)
  → Build chains (collapse non-branching same-power degree-2 edges)
  → Create half chains (forward/reverse half-edges)
  → Angular sort at vertices → sector coplanar pairs
  → Next pointers: connect half chains around each face
  → Find connected components → ray casting → seed exterior half chain (winding = 0)
  → BFS/DFS on half chain graph → winding number for every half chain
  → Filter by winding number
  → Dual cancellation + rebuild next pointers
  → Trace faces into rings → assemble multi-polygons
```

Complexity: O(n log n) for snap rounding (dominated by spatial index queries), O(n) for face traversal.

## Dependencies

- **Boost.Geometry**: rtree for spatial indexing, polygon validity, WKT I/O
- **Boost.Container**: flat_map for sorted associative containers
- **C++20**: `std::exclusive_scan`, abbreviated function templates

## Build

```bash
cmake -B build -G Ninja
cmake --build build
ctest --test-dir build
```

## References

1. Guibas, L. J., & Marimont, D. H. (1998). Rounding Arrangements Dynamically. *IJCGA*, 8(2), 157-176.
2. Greene, D. H., & Yao, F. F. (1986). Finite-resolution computational geometry. *FOCS 1986*.
3. Hobby, J. D. (1999). Practical segment intersection with finite precision output. *CGTA*, 13(4), 199-214.

---

## TODO

### Algorithm Completeness

- [ ] **Download Guibas & Marimont (1998) paper** and store in `docs/` directory for reference
- [ ] **Property merge**: merge polygon properties during boolean operations (like Boost.Polygon's property merge)

### Architecture

- [ ] **Graph-theoretic algorithms for planar graph operations**: All planar graph algorithms (connected components, face traversal, BFS propagation) should use standard graph theory algorithms, similar to Boost Graph Library style. All are O(n). Single-thread only — no multi-threading for now.

### Spatial Index

- [ ] **Replace Boost.Geometry rtree with custom uniform grid**: The spatial index for segment intersection and point-on-segment queries is the performance bottleneck. Implement a custom uniform grid spatial index instead of relying on Boost.Geometry rtree. This gives better cache locality and avoids Boost dependency for this part. Must include dedicated unit tests for the uniform grid.

### Performance Comparison

- [ ] **Add Clipper2 to CMake and benchmark head-to-head**: Integrate Clipper2 as a CMake dependency alongside the existing dependencies. Write benchmarks that run the same boolean operation test cases through both best_clipper and Clipper2, measuring and comparing runtime. The goal is to demonstrate that this algorithm outperforms Clipper2.
- [ ] **Benchmark comparison**: also compare against Boost.Polygon for performance
- [ ] **GPU acceleration**: use cuSpatial for spatial join (finding intersection points and segments on hot pixels), which is the most time-consuming step. Also explore migrating face traversal to GPU.
- [ ] **Bounding box grouping**: group segments by bounding box first, compute connected components before the full pipeline
- [ ] **Lower C++ standard requirement**: currently requires C++20; explore lowering to C++17 for wider compatibility

### Testing

- [ ] **Uniform grid spatial index unit tests**
- [ ] More comprehensive unit tests, especially for edge cases:
  - Coincident boundaries (holes exactly matching outer edges)
  - Self-intersecting input
  - Very large coordinates (overflow testing)
  - Degenerate inputs (zero-area polygons, collinear edges)
- [ ] Fuzz testing against Boost.Geometry's built-in boolean operations
