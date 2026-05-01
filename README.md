# best_clipper

A C++20 header-only library for polygon boolean operations (union, intersection, XOR, etc.) based on **snap rounding** and **winding number face traversal**.

## Core Ideas

### 1. Snap Rounding with Spatial Index

Snap rounding is a technique to convert arbitrary-precision geometric arrangements into integer-grid representations while preserving topological correctness. The theoretical foundation comes from:

> **Guibas, L. J., & Marimont, D. H. (1998).** *Rounding Arrangements Dynamically.* International Journal of Computational Geometry & Applications, Vol. 8, No. 2, pp. 157-176.

The key idea: all segment endpoints and intersection points are treated as **hot pixels** (integer grid points). Every segment is then "snapped" to pass through these hot pixels, producing a planar subdivision that is topologically equivalent to the true arrangement.

Unlike the original paper which handles dynamic (online) insertion, this library processes all segments in batch (offline). Instead of a dynamic rounding structure, a **spatial index tree** (Boost.Geometry rtree) is used to accelerate segment-intersection and point-on-segment queries. All geometric computations use **strict integer arithmetic** — no floating point inaccuracy.

Pipeline for snap rounding:
1. Find all segment-segment intersection points, snapped to integer grid (hot pixels)
2. For each original segment, find all hot pixels lying on it
3. Split each segment into sub-segments between consecutive hot pixels
4. Deduplicate edges (normalized so start < end)

### 2. Winding Number Boolean Operations

Boolean operations are defined strictly via **winding numbers** (绕数). For any point in the plane, its winding number counts the net number of counterclockwise boundary traversals enclosing it.

| Operation | Winding Number Condition |
|-----------|-------------------------|
| OR (union) | w > 0 |
| AND (intersection) | w > 1 |
| XOR | w = 1 |

#### Edge Power

Each edge has a **power**: the winding number of the face on its left side minus the winding number of the face on its right side. For an edge on a CCW outer boundary, power = +1; for CW (hole) boundary, power = -1. When edges from different input polygons overlap, their powers are **summed** (cancellation when equal and opposite). This is computed during snap rounding deduplication.

#### Chain Building

Adjacent edges with the **same power** and **no branching** (each intermediate vertex has exactly one incoming and one outgoing edge with the same power) are collapsed into a **chain**. This significantly reduces the graph size for subsequent steps.

#### Duplicated Chains (Half-Chain Model)

Each chain is split into two **directed half-chains** (duplicated chains, DCs), analogous to the half-edge data structure:
- **Forward DC** (even index): traverses the chain in its original direction; the face on the left of this DC is the "left face"
- **Reverse DC** (odd index): traverses the chain in reverse; the dual of the forward DC

The relationship: `dc.dual() = dc.id ^ 1` (XOR with 1 flips between forward and reverse).

The power of a directed chain equals the winding number difference when crossing from its left face to its right face. For a forward chain with power p, the forward DC has power p and the reverse DC has power -p.

#### Face Relationships Without Explicit Face IDs

Rather than computing explicit face IDs for each half-edge, the algorithm propagates winding numbers through three types of **coplanarity** (共面) relationships between DCs:

1. **Angular (sector) coplanarity**: At each vertex, all DCs starting from that vertex are sorted by CCW polar angle. Adjacent DCs in this sorted order share a face — the right face of one is the left face of the next. This gives: `dc_a.dual()` and `dc_b` are coplanar.

2. **Ray-casting coplanarity**: For each connected component, find the leftmost vertex. Cast a ray in the -x direction (leftward). At this vertex, the ray sits in the sector between two DCs in CCW order; the DC on the CW side of the ray is identified. Then find the nearest DC that this ray intersects. These two DCs are coplanar (same face). If the ray hits nothing, the leftmost DC's left face is the exterior (winding number = 0).

3. **Dual cancellation coplanarity**: When both dual DCs of a chain survive the winding number filter, they cancel each other and are removed. Their adjacent faces are then merged (the left face of one and the right face of the other become coplanar).

These coplanar pairs are fed into a **Disjoint Set Union (DSU)** to group DCs by face. Winding numbers are then propagated via BFS from known exterior faces (winding = 0). The key invariant: crossing a DC with power p changes the winding number by p (right face winding = left face winding + power).

#### Face Filtering and Output Construction

After computing winding numbers:
1. Filter DCs by winding number (e.g., w > 0 for union)
2. Perform dual cancellation: if both forward and reverse DCs of a chain survive, both are removed
3. Update "next" pointers: if a DC's next is deleted, follow `next[dual(deleted_dc)]` iteratively until a surviving DC is found (this guarantees termination because a surviving DC always exists in the cycle)
4. Trace surviving DCs along "next" pointers to form closed rings (polygons)
5. Classify rings as outer rings (CCW, positive area) or holes (CW, negative area) and assemble into multi-polygons

## Algorithm Pipeline Summary

```
Input polygons
  → Collect segments
  → Snap rounding (find intersections, split segments, deduplicate edges)
  → Assign edge powers (+1/-1)
  → Build chains (collapse non-branching degree-2 edges)
  → Create duplicated chains (forward/reverse half-chains)
  → Angular sort at vertices → sector coplanar pairs
  → Find connected components → ray casting → ray coplanar pairs
  → DSU + BFS propagation → winding numbers
  → Filter by winding number
  → Dual cancellation + update next pointers
  → Trace faces into rings → assemble multi-polygons
```

Complexity: O(n log n) for snap rounding (dominated by spatial index queries), O(n) for face traversal.

## Dependencies

- **Boost.Geometry**: rtree for spatial indexing, polygon validity, WKT I/O
- **Boost.Container**: flat_map for sorted associative containers
- **C++20**: ranges, concepts, `std::exclusive_scan`

## Build

```bash
cmake -B build -G Ninja
cmake --build build
ctest --test-dir build
```

## Why This Algorithm Is Theoretically Sound

This section explains why the algorithm has no inherent theoretical flaws — all potential issues are strictly implementation-level.

### Foundation 1: Snap Rounding Preserves Topology

The Guibas & Marimont (1998) paper proves that snapping all intersection points and endpoints to integer grid points (hot pixels) and re-routing segments through them produces a planar subdivision that is **topologically equivalent** to the true arrangement. Key guarantees:

- No two segments cross between grid lines without an intermediate hot pixel
- No segment passes too close to a hot pixel without passing through it
- The resulting planar graph is a valid subdivision of the plane

Since we batch process all segments (offline), the dynamic insertion complexity from the original paper is irrelevant. The spatial index tree (Boost.Geometry rtree) simply accelerates the geometric queries — it does not change the mathematical result.

### Foundation 2: Winding Numbers Are Mathematically Rigorous

For any point in a planar subdivision, its **winding number** with respect to a set of oriented closed curves is a well-defined integer. Within a single face of the subdivision, the winding number is **constant**. Therefore:

- **OR** = faces with w > 0 (point is inside at least one polygon)
- **AND** = faces with w > 1 (point is inside at least two polygons)
- **XOR** = faces with w = 1 (point is inside exactly one polygon)

This is the same mathematical principle behind the residue theorem in complex analysis. There is no ambiguity.

### Foundation 3: Edge Power Is Consistent and Additive

Each directed edge has a **power** = w(left face) - w(right face). For a CCW outer boundary edge, power = +1; for a CW hole boundary edge, power = -1.

When edges from different input polygons overlap (coincident boundaries), their powers **sum**. If the sum is zero, the edge carries no net winding contribution and is **safely removed** — both adjacent faces have the same winding number, and removing the edge does not change any topological property of the subdivision.

### Foundation 4: Chain Building Is Lossless

Collapsing a sequence of degree-2 non-branching edges with identical power into a single chain is an information-preserving compression. Each vertex along the chain has exactly one incoming and one outgoing edge with the same power, so no topological ambiguity is introduced. The chain endpoints are exactly the branch points of the graph.

### Foundation 5: Three Sources of Coplanarity Are Necessary and Sufficient

The face structure of the planar subdivision is fully determined by three types of coplanarity (共面) relations between directed chains (DCs):

1. **Angular (sector) coplanarity**: At each vertex, DCs are sorted by CCW polar angle. In a planar subdivision, the right face of one departing DC IS the left face of the next departing DC in CCW order. Therefore `dc_a.dual()` and `dc_b` (adjacent in sorted order) necessarily belong to the same face. This is a geometric fact, not a heuristic.

2. **Ray-casting coplanarity**: From the leftmost vertex of each connected component, a -x ray identifies the face on the exterior side of the component. The DC on the CW side of the ray and the nearest DC intersected by the ray belong to the same face. If no intersection exists, the leftmost DC's left face is the exterior (w = 0). This determines the **reference frame** for winding numbers.

3. **Dual cancellation coplanarity**: When both forward and reverse DCs of a chain survive filtering, they form an interior boundary that separates two faces with the same winding number. Removing both and merging their adjacent faces is correct because the chain's own power cancels out (p + (-p) = 0).

These three sources cover **all** coplanarity relations. Together with the fact that every edge connects exactly two faces, the face connectivity graph is fully determined.

### Foundation 6: Winding Number Propagation Is Complete

Starting from the exterior face (w = 0, identified by ray casting), BFS along face-adjacency edges propagates winding numbers. The invariant: crossing a DC with power p changes the winding number by p. Since the face adjacency graph is connected (any face is reachable from the exterior through a finite sequence of boundary crossings), every face receives a correct winding number. There are no "unreachable" faces in a correct planar subdivision.

### Foundation 7: Dual Cancellation "Next" Update Terminates

When a DC is deleted (dual cancellation), its "next" pointer must be updated. For each DC `p` with `next[p] == dead_dc`, the new next is `next[dual(dead_dc)]`. If that next is also dead, continue following `next[dual(dead_next)]` iteratively. This chain of indirection always terminates at a surviving DC because: (a) the graph is finite, (b) the cycle of DCs around each vertex never becomes empty (otherwise the vertex would be isolated and should not have survived).

### Foundation 8: Self-Intersecting Input Is Handled

Self-intersecting input polygons are handled naturally: the snap rounding step splits all intersecting segments at their intersection points, creating a proper planar subdivision. The winding number of each resulting face correctly reflects the net enclosure count. No special pre-processing beyond `bg::correct()` is needed.

### Summary

Every step of the algorithm is a direct consequence of planar subdivision topology. There are no approximations, no heuristics, and no special cases that require exceptions. The algorithm is correct **in theory**.

Any bugs that exist are **implementation errors** — code that does not faithfully execute the mathematical steps described above. The known implementation issues are listed in the TODO section below.

## References

1. Guibas, L. J., & Marimont, D. H. (1998). Rounding Arrangements Dynamically. *IJCGA*, 8(2), 157-176.
2. Greene, D. H., & Yao, F. F. (1986). Finite-resolution computational geometry. *FOCS 1986*.
3. Hobby, J. D. (1999). Practical segment intersection with finite precision output. *CGTA*, 13(4), 199-214.

---

## TODO

### Algorithm Completeness

- [ ] **Download Guibas & Marimont (1998) paper** and store in `docs/` directory for reference
- [ ] **XOR operation**: currently only union (w > 0), intersection (w > 1), and self-union are implemented. XOR (w = 1) needs to be added
- [ ] **Difference operation**: A - B (faces where w_A > 0 and w_B = 0) is not yet implemented
- [ ] **Property merge**: merge polygon properties during boolean operations (like Boost.Polygon's property merge)

### "Next" Pointer Update After Dual Cancellation

- [ ] The current implementation rebuilds next/prev pointers from scratch by re-scanning sorted_dcs at each vertex after dual cancellation. This is incorrect. The correct approach: when a DC `d` is deleted (dead), for each DC `p` whose `next[p] == d`, update `next[p] = next[dual(d)]`. If that next is also dead, continue following `next[dual(dead_dc)]` iteratively until a surviving DC is found. This chain of indirection is guaranteed to terminate because at least one surviving DC exists in the cycle. The same logic applies to prev pointers.

### Implementation Issues (Not Theoretical Bugs)

- [ ] **Zero-power edges should be removed**: `std::erase_if` for zero-power edges is commented out (line 268). Coincident boundaries cancel via power summation (winding number semantics), so zero-power edges carry no net winding contribution and should be deleted. They should not be kept.
- [ ] **Unreachable-face fallback should not exist** (line 686-687): `if (q.empty()) { rwind[root[0]] = 0; ... }` — all faces should be reachable from exterior faces via BFS propagation. If the queue is empty, there is a bug in coplanarity propagation, and silently defaulting to winding 0 masks the real problem. This fallback should be removed or replaced with an assertion/error.

### Architecture

- [ ] **Graph-theoretic algorithms for planar graph operations**: All planar graph algorithms (connected components, face traversal, BFS propagation) should use standard graph theory algorithms, similar to Boost Graph Library style. All are O(n). Single-thread only — no multi-threading for now.

### Spatial Index

- [ ] **Replace Boost.Geometry rtree with custom uniform grid**: The spatial index for segment intersection and point-on-segment queries is the performance bottleneck. Implement a custom uniform grid spatial index instead of relying on Boost.Geometry rtree. This gives better cache locality and avoids Boost dependency for this part. Must include dedicated unit tests for the uniform grid.

### Performance Comparison

- [ ] **Add Clipper2 to CMake and benchmark head-to-head**: Integrate Clipper2 as a CMake dependency alongside the existing dependencies. Write benchmarks that run the same boolean operation test cases through both best_clipper and Clipper2, measuring and comparing runtime. The goal is to demonstrate that this algorithm outperforms Clipper2.
- [ ] **Benchmark comparison**: also compare against Boost.Polygon for performance
- [ ] **GPU acceleration**: use cuSpatial for spatial join (finding intersection points and segments on hot pixels), which is the most time-consuming step. Also explore migrating face traversal to GPU.
- [ ] **Multi-precision integer support**: use GMP or similar for exact arithmetic in snap rounding (currently using `int` with overflow potential for large coordinates)
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
