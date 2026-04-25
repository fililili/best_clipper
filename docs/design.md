# best_clipper Design Document

## Overview

`best_clipper` is a C++20 header-only library for polygon boolean operations (union, intersection) using snap rounding planar graph construction followed by chain-based face traversal. The algorithm leverages the theoretical foundation of snap rounding to guarantee topological correctness while maintaining O(n log n) time complexity.

## Theoretical Foundation

### Snap Rounding

The algorithm follows the snap rounding framework introduced by Greene & Hobby and rigorously analyzed in:

> **Guibas, L. J., & Marimont, D. H. (1998).** *Rounding Arrangements Dynamically.* International Journal of Computational Geometry & Applications, Vol. 8, No. 2, pp. 157–176. DOI: [10.1142/S0218195998000096](https://doi.org/10.1142/S0218195998000096)

The paper proves that snapping intersection points to integer grid points (hot pixels) preserves the topological structure of the arrangement. Key lemma: if all hot pixels are added as vertices and segments are re-routed through them, the resulting planar subdivision is topologically equivalent to the true arrangement.

### Winding Number

The winding number (绕数) of a face is computed using the residue theorem. For a point in a face, the winding number counts the net number of counterclockwise turns around the point by the boundary. For boolean operations:

- **Union**: faces with winding number > 0 survive (point is inside at least one polygon)
- **Intersection**: faces with winding number > 1 survive (point is inside both polygons)

### Chain-based Face Traversal

Rather than operating on individual half-edges, we collapse non-branching edges into **chains** to reduce the graph size. This is the key innovation over standard half-edge approaches.

## Algorithm Pipeline (14 Steps)

### Step 1: Find All Intersection Points (Hot Pixels)

Input segments are decomposed from multi-polygons. All pairwise segment intersections are computed. Each intersection point is snapped to the nearest integer grid point (hot pixel). All segment endpoints are also added as hot pixels.

### Step 2: Split Segments at Hot Pixels

For each segment, find all hot pixels that lie on it. Sort these pixels along the segment direction. Form new sub-segments between consecutive hot pixels.

### Step 3: Deduplicate Edges and Compute Power

Each edge is normalized so that `start < end`. Duplicate edges (same start/end) have their power summed. The power of an edge is +1 for a forward traversal (from the left polygon boundary) and -1 for backward. The power represents the winding number contribution when crossing the edge from left to right.

Edges with power == 0 (equal forward and backward traversals) are removed.

### Step 4: Build Chains

Consecutive edges where each vertex has exactly one incoming and one outgoing edge with the same power are merged into chains. Chains reduce the number of elements for subsequent steps.

A vertex is a chain endpoint if:
- Its out_degree != 1, or
- Its in_degree != 1, or  
- out_power != in_power

Starting from each endpoint, follow edges until another endpoint is reached.

Pure cycles (where all vertices satisfy the chain-continuation condition) are handled separately.

### Step 5: Create Duplicated Chains

Each chain produces two directed half-chains (dual chains / duplicated chains):
- **Forward** (even index): traverses the chain in original direction, belongs to the left face
- **Reverse** (odd index): traverses the chain in reverse direction, belongs to the right face

The power of a duplicated chain is the winding number difference when crossing from its left face to its right face.

### Step 6: Angular Sort at Vertices

At each vertex, duplicated chains starting from that vertex are sorted counterclockwise by their direction vector. Adjacent duplicated chains in this sorted order share a face — the face on the right side of one is the face on the left side of the next.

### Step 7: Face Co-planarity Propagation

Given the angular sort, adjacent duplicated chains a and b establish that a's dual and b are coplanar (belong to the same face) with the same winding number.

### Step 8: Connected Components

Vertices are grouped into connected components based on shared duplicated chains. Within each component, face winding numbers are related.

### Step 9: Ray Casting for Exterior Face

For each connected component, find the leftmost vertex. Cast a ray in the -x direction. Find the nearest intersecting duplicated chain. The face on the left side of this chain belongs to the same connected component. If the ray hits no chain, the leftmost vertex's left face is the exterior (winding number = 0).

### Step 10: Compute All Winding Numbers

Using DFS/BFS from the known exterior face(s), propagate winding numbers through all coplanar relationships. The difference in winding number between adjacent faces equals the power of the separating duplicated chain.

### Step 11: Filter by Winding Number

Duplicated chains are filtered based on their face's winding number:
- Union: keep chains whose left face has winding number > 0
- Intersection: keep chains whose left face has winding number > 1

### Step 12: Dual Cancellation

If both directions of a chain survive, they cancel out (remove from output). Cancellation also modifies the adjacency of neighboring duplicated chains: the predecessor of one is linked to the successor of the other, and vice versa.

### Step 13: Form Rings

Surviving duplicated chains are connected end-to-end using the angular ordering at each vertex. Each connected cycle forms a ring.

### Step 14: Assemble Multi-Polygons

Rings are classified as:
- **Outer rings** (CW, positive area): the ring with minimum x-coordinate is the outer ring
- **Holes** (CCW, negative area): nested inside outer rings

Holes are assigned to their parent outer ring using point-in-polygon tests (DE-9IM `TFF` mask via Boost.Geometry rtree).

## Dependencies

- **Boost.Geometry**: rtree for spatial indexing, polygon validity, WKT I/O
- **Boost.Container**: flat_map for sorted associative containers
- **C++20**: ranges, concepts, std::exclusive_scan

## References

1. Guibas, L. J., & Marimont, D. H. (1998). Rounding Arrangements Dynamically. *IJCGA*, 8(2), 157-176.
2. Greene, D. H., & Yao, F. F. (1986). Finite-resolution computational geometry. *FOCS 1986*.
3. Hobby, J. D. (1999). Practical segment intersection with finite precision output. *CGTA*, 13(4), 199-214.
4. CGAL 2D Snap Rounding documentation: https://doc.cgal.org/latest/Snap_rounding_2/
