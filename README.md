# best_clipper
A new algorithm to do polygon boolean operation like the library clipper2. In the future, may support the property merge function (boost polygon has implementation)
The core idea of this library is use planar face traversal algorithm, which is an O(n) function and easy to understand. So the code is very simple and fast. And to build a planar graph, I will use snap rounding algorithm, snap rounding will use exsiting spatial index tree implementaion (boost geometry rtree), the time complexity is nearly O(nlogn). (n is the number of interected points and endpoints)
Shortly speaking, the algorithm has two step:
1. do snap rounding to bulid a planar graph.O(nlogn)
2. do planar face traversal. O(n)

I have not found similar algorithm before, if anybody find one, please cantact me.

A basic code and test have been write already, still need do:

0. not benchmark compare now, at current testing, the implementation is fast (3s for 1 000 000 vertices in myself computer.) But have not compared with other implementation now.
1. after merege faces, they may co-point, need to find a goode algorithm, and the find parent hole progress maybe time cost, need find a simple way.
3. try to use cuSpatial to do spatial join (find intersections points and find segs on pixels, which cost most of the time), also, find the way to migrate other part of algorithms to GPU, like planer face traversal.
4. use a precise big number library (like gmp) to write the snap rounding realting geometry algorithms.
5. add more unit test, and migrate to gtest framework.
6. it will be a head only lib based on other libs, but will try use cmake to organize the lib.
7. may group by bounding box first and calculate graph connected component firstly.
8. only support add, intersection now, need support other boolean operation and test.
9. down the c++ standard.
