# best_clipper

A C++20 library for polygon boolean operations (union, intersection, XOR, etc.) based on **snap rounding** and **winding number face traversal**.

## Core Ideas

### 1. Glocal Snap Rounding with Spatial Index

Global snap rounding is a technique to convert arbitrary-precision geometric arrangements into integer-grid representations while preserving topological correctness. Global means it will take **every** segments into consideration when snaping every**one** point. And the snap will not cross any segment, or call it topo inconsistency. The theoretical foundation comes from:

> **Guibas, L. J., & Marimont, D. H. (1998).** *Rounding Arrangements Dynamically.* International Journal of Computational Geometry & Applications, Vol. 8, No. 2, pp. 157-176.

The key idea: all segment endpoints and intersection points are treated as **hot pixels** (integer grid points). Every segment is then "snapped" to passed hot pixels.

Unlike the original paper which handles dynamic (online) insertion, this library processes all segments in batch (offline). Instead of a dynamic rounding structure, a **spatial index tree** (uniform grid now) is used to accelerate segment-intersection and point-on-segment queries. All geometric computations use **strict integer arithmetic** — no floating point inaccuracy.

Pipeline for snap rounding:
1. Find all segment-segment intersection points, snapped to integer as hot pixels
2. All endpoint of segment is hot pixels, too.
3. sort and unique hot pixels.
4. For each original segment, find all hot pixels lying on it
5. Split each segment into sub-segments between consecutive hot pixels
6. merge sub-segments and build chains.

### 2. Winding Number

**Definition.** For a point P in the plane, its **winding number** w(P) is defined by using **Residue theorem**. According to **Residue theorem**, even point in same face have same **winding number**. PS: A face is a maximal connected region.
### Edge Power

An edge split two face, so we can view an edge as two **half edge**, **Forward half edge** and **Reverse half edge**. Each half edge belong to a face. Normally, people call it **half-edge data structure**. What's more, we call **Forward half edge** and **Reverse half edge** as **dual half edge**.
Then we can define edge power as:
**winding number** of the face that **Forward half edge** belong to
sub
**winding number** of the face that **Reverse half edge** belong to
= power of edge
According to the definition, we can simplify merge edge power by add.

### Chain Building

Adjacent edges with the **same power** and **no branching** (each intermediate vertex has exactly one incoming and one outgoing edge with the same power) are collapsed into a **chain**. This significantly reduces the graph size for subsequent steps.

### Half Chains (Half-Edge Data Structure)

Same as edge, each chain is split into two **half chains**, **Forward half chain** and **Reverse half chain**. :

### Following input orientation
In above definition, we only said: **Forward half edge** and **Reverse half edge**. Each half edge belong to a face.
However, we don't clearly define which face **Forward half edge** or **Reverse half edge** belong to. The definition will following a rule:
For a valid multipolygon, even interior point's winding number is 1. Therefore, **Forward half edge** belong to interior face, and **Reverse half edge** belonging to outer face.

| Operation | Winding Number Condition |
|-----------|-------------------------|
| OR (union) | w > 0 |
| AND (intersection) | w > 1 |
| XOR | w == 1 |
| A - B | w = w_A + (-1) * w_B > 0 |

At current, we use **Boost.Geometry multipolygon set**. For boost, outer ring match CW, so every half-edge belongs to the face on its right. In the future, we will support CCW, too. For multipolygon that has CCW outer ring, every half-edge belongs to the face on its left. (we said left and right for half edge that is from down to up.)
A face is contained by many half chains. And we can find the **next half chains** by sorting half chains that has same start point.
However, same face have many rings, we can find them by ray cast. **Ray cast** can also check whether a half chain is belong to the exterior face whose winding number is always zero.
By **Dual half chain**, **next half chain**, **Ray cast pair half chain** and **Ray cast exterior half chain**, we can use DFS to calculte the winding number of every half chain and then filter them.

### Output Construction

1. Filter half chains: keep those whose `winding` satisfies the operation (e.g. `winding > 0` for union)
2. Dual cancellation: if both forward and reverse half chains of a chain survive, kill both and merge faces.
3. Rebuild `next_half_chain` respecting surviving half chains and face boundaries
4. Trace surviving half chains along `next_half_chain` into closed rings.
5. Assemble rings into polygons. **next half chain**, **Ray cast pair half chain** and **merged dual half** belong to same face, a face is a polygon.
6. left most point belong to outer ring, rest = holes

#### Self-Intersecting Input

Self-intersecting input polygons are handled naturally: the snap rounding step splits all intersecting edges at their intersection points, creating a proper planar subdivision.

## Complexity
O(n log n) for snap rounding (dominated by spatial index queries), O(n) for face traversal.

## Dependencies

- **Boost.Geometry**: rtree for spatial indexing, polygon validity, WKT I/O
- **clipper2**: for testing.

## Build

```bash
cd best_cliiper && mkdir build && cmake -B build
```

## References

1. Guibas, L. J., & Marimont, D. H. (1998). Rounding Arrangements Dynamically. *IJCGA*, 8(2), 157-176.
2. Greene, D. H., & Yao, F. F. (1986). Finite-resolution computational geometry. *FOCS 1986*.
3. Hobby, J. D. (1999). Practical segment intersection with finite precision output. *CGTA*, 13(4), 199-214.

## TODO
1. support two stage cluster, for stage 1, use bbox. And totally sperate different cluster. For stage 2, use bbox and segs, different cluster has different snap rounding and chains build logic, but share ray casting logic.
2. support connected_component with rank.
3. For two dual half edge of an edge, only one half edge is possible for ray casting, so ray cast can store half edge index directly.
4. use power_2 to handle robust input.
5. sizing support
6. clipping linestring support (need study to define new power)
7. figure out interger overflow risk
8. add more practice testing
9. sort ray query's for cache friendly

## More study
1. CUDA speed up for both graph and BVH glgo.
2. More quickly snap rounding geometry functions.
3. Multi thread support for both graph and BVH algo.

### Performance Comparison
Claude code auto generated code, need check and add more practice cases.

### Testing
Claude code auto generated code, need check and add more.
