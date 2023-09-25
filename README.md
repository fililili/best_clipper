# best_clipper
A new algorithm to do polygon boolean operation like the library clipper2. In the future, may support the property merge function (boost polygon has implementation)
The core idea of this library is use planar face traversal algorithm, which is an O(n) function and easy to understand. So the code is very simple and fast. And to build a planar graph, I will use snap rounding algorithm, snap rounding will use exsiting spatial index tree implementaion (boost geometry rtree), the time complexity is nearly O(nlogn).
Shortly speaking, the algorithm has two step:
1. do snap rounding to bulid a planar graph.
2. do planar face traversal.
I have not found similar algorithm before, if anybody find one, please cantact me.

A basic code and test have been write, still need do:
0. not benchmark compare now, at current testing, the implementation is fast (10s for 1 000 000 vertices in myself computer.) But have not compared with other implementation now.
1. the planar face traversal implementaion is a little complex, so bug may exist. Need to find a way to simple the graph data structure.
2. have not handle co-point situation, and the find parent hole progress may be time cost, need find a simple way.
3. try to use cuSpatial to do spatial join (find intersections points and find segs on pixels, which cost most of the time), also, find the way to migrate other part of algorithms to GPU, like planer face traversal.
4. use a precise big number library (like gmp) to write the snap rounding realting geometry algorithms.
5. add more unit test, and migrate to gtest framework.
6. it will be a head only lib based on other libs, but will try use cmake to organize the lib.
