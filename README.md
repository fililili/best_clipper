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

#### Half Chains (Half-Edge Data Structure)

Each chain is split into two **half chains** (HCs), which form a standard **half-edge data structure**:

- **Forward HC** (even index): traverses the chain in its original direction, belonging to the face on its **left** side
- **Reverse HC** (odd index): traverses the chain in the opposite direction, belonging to the face on its **left** side (the other face of the chain)

This is the standard half-edge convention: every half-edge belongs to the face on its left. The key operations:

| Half-Edge Concept | Implementation |
|------------------|----------------|
| twin / opposite | `hc.dual()` = `hc.id ^ 1` |
| next (around face) | `next_hc[hc.id]` = next half-edge in the same face |
| power (winding delta) | forward: `+chain_power`, reverse: `-chain_power` |

#### Face Determination via Coplanarity

Faces are not tracked explicitly. Instead, half chains are grouped into faces through three **coplanarity** (共面) relations:

1. **Sector coplanarity**: At each vertex, departing HCs are sorted by CCW angle. Between adjacent HCs `prev` and `cur`, the sector lies on the left of `prev` and on the right of `cur`. Since a half-edge belongs to its left face, `prev` and `cur.dual()` belong to the same face. Hence `coplanar(prev, cur.dual())`.

2. **Ray-casting coplanarity**: From the leftmost vertex of each component, a -x ray determines the exterior face. The HC whose left face contains the ray and the nearest HC hit by the ray are coplanar. If the ray hits nothing, that HC's left face is the exterior (winding = 0).

3. **Dual cancellation coplanarity**: When both forward and reverse HCs of a chain survive filtering, they cancel each other. Their adjacent faces merge: `coplanar(fwd, rev)`.

Coplanar pairs feed a **Disjoint Set Union (DSU)**, grouping HCs by face. Next pointers link HCs around each face boundary. The face adjacency of the planar subdivision is then: crossing a HC with power p changes winding by p.

#### Face Filtering and Output Construction

After computing winding numbers:
1. Filter HCs by winding number (e.g., w > 0 for union) → surviving HCs
2. Dual cancellation: if both forward and reverse HCs of a chain survive, both are removed and their adjacent faces are merged
3. Rebuild `next_hc` pointers respecting new face boundaries
4. Trace surviving HCs along `next_hc` to form closed rings
5. Assemble rings into polygons (outer ring by max area, remaining become holes), orienting outer CCW and holes CW as required by Boost.Geometry

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
  → Find connected components → ray casting → ray coplanar pairs
  → DSU + BFS propagation → winding numbers
  → Filter by winding number
  → Dual cancellation + rebuild next pointers
  → Trace faces into rings → orient → assemble multi-polygons
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

The face structure of the planar subdivision is fully determined by three types of coplanarity (共面) relations between half chains:

1. **Sector coplanarity**: At each vertex, departing HCs are sorted by CCW angle. The sector between adjacent HCs `prev` and `cur` lies on the left of `prev` and on the right of `cur`. Since a half-edge belongs to its left face, `prev` and `cur.dual()` belong to the same face. This is a geometric fact, not a heuristic.

2. **Ray-casting coplanarity**: From the leftmost vertex of each connected component, a -x ray identifies the face on the exterior side. The HC whose left face contains the ray and the nearest HC hit by the ray are coplanar. If no hit, that HC's left face is the exterior (w = 0). This establishes the **reference frame** for winding numbers.

3. **Dual cancellation coplanarity**: When both forward and reverse HCs of a chain survive filtering, they form an interior boundary separating two faces with the same winding number. Merging these faces is correct because the chain's power cancels out (p + (-p) = 0).

These three sources cover **all** coplanarity relations in the planar subdivision.

### Foundation 6: Winding Number Propagation Is Complete

Starting from the exterior face (w = 0, identified by ray casting), BFS along face-adjacency edges propagates winding numbers. The invariant: crossing a HC with power p changes the winding number by p. Since the face adjacency graph is connected (any face is reachable from the exterior through a finite sequence of boundary crossings), every face receives a correct winding number. There are no "unreachable" faces in a correct planar subdivision.

### Foundation 7: Dual Cancellation "Next" Update Terminates

When both forward and reverse HCs of a chain are removed (dual cancellation), their `next_hc` pointers are bypassed. For each HC `p` whose `next_hc[p]` points to a dead HC `d`, the new next is `next_hc[d.dual()]`. If that is also dead, the chain of indirection is followed until a surviving HC is reached. This always terminates because the cycle around each vertex never becomes empty.

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

- [ ] The current implementation rebuilds next/prev pointers from scratch by re-scanning sorted_hcs at each vertex after dual cancellation. This is incorrect. The correct approach: when a HC `d` is deleted (dead), for each HC `p` whose `next[p] == d`, update `next[p] = next[dual(d)]`. If that next is also dead, continue following `next[dual(dead_dc)]` iteratively until a surviving HC is found. This chain of indirection is guaranteed to terminate because at least one surviving HC exists in the cycle. The same logic applies to prev pointers.

### Implementation Issues (Not Theoretical Bugs)

- [x] ~~**Zero-power edges should be removed**~~ ✓ Fixed
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
