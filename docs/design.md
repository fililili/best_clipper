# best_clipper Design Document

## Overview

`best_clipper` is a C++20 header-only library for polygon boolean operations (union, intersection, XOR, difference, self-union) using snap rounding planar graph construction followed by winding number face traversal. The algorithm leverages the theoretical foundation of snap rounding to guarantee topological correctness while maintaining O(n log n) time complexity.

## Theoretical Foundation

### Snap Rounding

The algorithm follows the snap rounding framework introduced by Greene & Hobby and rigorously analyzed in:

> **Guibas, L. J., & Marimont, D. H. (1998).** *Rounding Arrangements Dynamically.* International Journal of Computational Geometry & Applications, Vol. 8, No. 2, pp. 157–176. DOI: [10.1142/S0218195998000096](https://doi.org/10.1142/S0218195998000096)

The paper proves that snapping intersection points to integer grid points (hot pixels) preserves the topological structure of the arrangement. Key lemma: if all hot pixels are added as vertices and segments are re-routed through them, the resulting planar subdivision is topologically equivalent to the true arrangement.

### Winding Number

The winding number (绕数) of a face is computed using the residue theorem. For a point in a face, the winding number counts the net number of counterclockwise turns around the point by the boundary. For boolean operations:

- **Union**: faces with winding number > 0 survive
- **Intersection**: faces with winding number > 1 survive
- **XOR**: faces with winding number == 1 survive
- **Difference (A - B)**: reverse B's segments, then faces with winding number > 0 survive
- **Self-union** (`self_or`): union of a polygon with itself, resolves self-intersections

### Chain-based Face Traversal

Rather than operating on individual half-edges, we collapse non-branching edges into **chains** to reduce the graph size. This is the key innovation over standard half-edge approaches.

## Algorithm Pipeline

### Steps 1-3: Edges with Power

All segment endpoints and pairwise segment intersections are snapped to integer grid points (**hot pixels**). Each segment is split into sub-segments between consecutive hot pixels lying on it. Edges are normalized (`start < end`) and deduplicated: overlapping edges from different input polygons have their powers summed. Power = +1 when an edge goes from lower to higher pixel index, -1 when reversed. Zero-power edges (equal forward and reverse traversal) are removed.

### Step 4: Build Chains

Consecutive edges where each intermediate vertex has exactly one incoming and one outgoing edge with the same power are merged into **chains**. This is a lossless compression: chain endpoints are exactly the branch points of the graph, and all topological information is preserved.

### Steps 5-6: Half-Chain Graph

Each chain produces two **half chains** following Boost.Geometry's right-face convention:
- **Forward** (even index): belongs to the face on its **right** side
- **Reverse** (odd index): belongs to the face on its **right** side (the other face)

At each vertex, departing half chains are sorted by CCW angle. The sector between adjacent half chains `prev` and `cur` (CCW) lies on the right of `cur`. Since a half chain belongs to its right face, `cur` and `prev.dual()` are **coplanar** (same face). Next pointers connect half chains around each face: `next[prev.dual()] = cur`.

### Step 7: Ray Casting

For each connected component, the leftmost vertex is identified. A ray is cast in the -x direction. The half chain whose right face contains the ray, and the nearest half chain hit by the ray, are coplanar. If nothing is hit, that half chain's right face is the exterior (winding = 0). This seeds winding number propagation.

### Step 8: Compute Winding Numbers

A weighted graph on half-chain IDs is built: coplanar/ray pairs (weight 0) and dual edges `i ↔ i^1` (weight `-power(i)`). A DFS from the exterior seeds (w = 0) propagates winding numbers to all reachable half chains.

### Steps 9-10: Build Output Polygons

1. Filter half chains by winding number according to the operation
2. Dual cancellation: if both forward and reverse half chains of a chain survive, kill both
3. Rebuild next pointers, bypassing dead half chains via indirection
4. Trace surviving half chains along next pointers into closed rings
5. Assemble rings into polygons (ring with minimum x-coordinate = outer, rest = holes)

## Dependencies

- **Boost.Geometry**: rtree for spatial indexing, segment/polygon model types
- **Boost.Container**: flat_map for sorted associative containers
- **C++20**: `std::exclusive_scan`, abbreviated function templates

## References

1. Guibas, L. J., & Marimont, D. H. (1998). Rounding Arrangements Dynamically. *IJCGA*, 8(2), 157-176.
2. Greene, D. H., & Yao, F. F. (1986). Finite-resolution computational geometry. *FOCS 1986*.
3. Hobby, J. D. (1999). Practical segment intersection with finite precision output. *CGTA*, 13(4), 199-214.
4. CGAL 2D Snap Rounding documentation: https://doc.cgal.org/latest/Snap_rounding_2/
