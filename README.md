# best_clipper
A new algorithm to do polygon boolean operation like the library clipper2. In the future, may support the property merge function (boost polygon has implementation)
The core idea of this library is use planar face traversal algorithm, which is an O(n) function and easy to understand. So the code is very simple and fast. I have not found similar algorithm before, if anybody find one, please cantact me.
Simplify speaking, will do snap rounding firstly, then do planar face traversal.
A basic code and test have been write, still need do:
1. bug may exist, need exatrct some code to handle my special flat map and faces by directed edge.
2. try to use cuSpatial to do spatial join, also, find the way to migrate other part of algorithms to GPU.
3. use a precise big number library (like gmp) to write the snap rounding realting geometry algorithms.
4. add more unit test, and migrate to gtest framework.
5. it will be a head only lib based on other libs, but will try use cmake to organize the lib.
