// Copyright 2018 Toyota Research Institute.  All rights reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

template <typename T> using VecXt = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using MatXt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// Defines a collection of linear algebra utils.
namespace linalg_utils {

// Generates `N` linearly spaced values between `lb` and `ub`, inclusive.
template <typename T>
std::vector<T> Linspace(const T lb, const T ub, const int N) {
  if (lb == ub && N == 1) {
    std::vector<T> vec;
    vec.push_back(lb);
    return vec;
  } else if (N <= 1)
  {
    std::cout << "N: " << N << " upper bound: " << ub << " lower bound: " << lb << std::endl;
    throw std::invalid_argument("[linalg_utils::Linspace] Number of discretization must be greater than 1 for unequal lower and upper bounds!");
  }
  if (lb >= ub)
    throw std::invalid_argument("[linalg_utils::Linspace] lower bound must be larger (or equal to) higher bound!");

  // Generate the evenly spaced vector.
  T h = (ub - lb) / static_cast<T>(N - 1);
  std::vector<T> vec;
  for (int n = 0; n < N; n++) {
    vec.push_back(lb + n*h);
  }
  return vec;
}

// Generates a linearly spaced vector between `lb` and `ub` with spacing of `h`, inclusive.
template <typename T>
std::vector<T> Arange(const T lb, const T ub, const float h = 1.0, bool endpoint = true)
{
  if (ub - lb < h)
    throw std::invalid_argument("[linalg_utils::Arange] discretization is larger than the provided range");

  std::vector<T> vec;

  // Do this instead of adding to avoid rounding errors
  int N = std::round((ub - lb) / h);

  // Add the last value if desired, only if it is not already added by the loop. Should return N+1 if this is enabled.
  if (endpoint) N+=1;

  for (int n = 0; n < N; n++) {
    vec.push_back(lb + n*h);
  }

  return vec;
}

// Generates a linearly spaced vector between `lb` and `ub` with spacing of `h`, inclusive. If ub < lb, returns an empty vector. Else, returns.
inline Eigen::VectorXd ArangeEigenUnsigned(const double lb, const double ub, const double h = 1.0, bool endpoint = true)
{
  if (ub - lb < h)
    throw std::invalid_argument("[linalg_utils::Arange] discretization is larger than the provided range");
  if (ub < lb) {
    Eigen::VectorXd vec(0);
    return vec;
  }
  if (ub == lb) {
    Eigen::VectorXd vec(1);
    vec(0) = ub;
    return vec;
  }

  // Do this instead of adding to avoid rounding errors
  int N = std::round((ub - lb) / h);


  // Add the last value if desired, only if it is not already added by the loop. Should return N+1 if this is enabled.
  if (endpoint) N+=1;

  Eigen::VectorXd vec(N);
  for (int n = 0; n < N; n++) {
    vec(n) = (lb + n*h);
  }

  return vec;
}

inline Eigen::VectorXi ArangeEigenUnsigned(const int lb, const int ub, const int h = 1, bool endpoint = true)
{
  if (ub < lb) {
    Eigen::VectorXi vec(0);
    return vec;
  }
  if (ub == lb) {
    Eigen::VectorXi vec(1);
    vec(0) = ub;
    return vec;
  }

  // Do this instead of adding to avoid rounding errors
  int N = std::round((ub - lb) / h);


  // Add the last value if desired, only if it is not already added by the loop. Should return N+1 if this is enabled.
  if (endpoint) N+=1;

  Eigen::VectorXi vec(N);
  for (int n = 0; n < N; n++) {
    vec(n) = (std::round(lb + n*h));
  }

  return vec;
}

// Cumtrapz computes the approximate cumulative integral of Y via the
// trapezoidal method with unit spacing. Analogous to MATLAB's cumtrapz.
template <typename T>
VecXt<T> Cumtrapz(T dx, const VecXt<T>& f)
{

  VecXt<T> y;
  y.setZero(f.size(), 1);

  y(0) = dx * f(0);
  for (size_t i = 1; i < y.size(); ++i) {
      y(i) = dx * f(i) + y(i-1);
  }

  return y;
}

// Cumtrapz computes the approximate cumulative integral of Y via the
// trapezoidal method with unit spacing. Analogous to MATLAB's cumtrapz.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
Container<T> Cumtrapz(T dx, const Container<T>& f)
{

  Container<T> y;
  y.resize(f.size(), 1);

  y[0] = dx * f[0];
  for (size_t i = 1; i < y.size(); ++i) {
      y[i] = dx * f[i] + y[i-1];
  }

  return y;
}

// Normalizes a heading that is between 0 and 2 pi
// into a heading that is between -pi and pi.
template<typename T>
T NormalizeHeading(T a)
{
  while (a > M_PI) {
    a -= 2 * M_PI;
  }

  while (a < -M_PI) {
    a += 2 * M_PI;
  }

  return a;
}

// Computes the radian difference in heading
// (a - b). The result will be between -pi and pi.
// Pre conditions: a and b are both between -pi and pi.
template<typename T>
T HeadingDifference(T a, T b)
{
  T diff = a - b;

  if (diff > M_PI)
  {
    return diff - 2 * M_PI;
  }
  else if (diff < -M_PI)
  {
    return diff + 2 * M_PI;
  }

  return diff;
}

// Closest point on a line:
// project vector AP onto vector AB, then add the resulting vector to point A
// A + dot(AP,AB) / dot(AB,AB) * AB
inline Eigen::Vector3d ClosestPointOnALine(const Eigen::Vector3d& origin, const Eigen::Vector3d& point, const Eigen::Vector3d& line_endpoint)
{
  Eigen::Vector3d vector = point - origin;
  Eigen::Vector3d line = line_endpoint - origin;
  double scale = vector.dot(line)/line.dot(line);
  return origin + line * scale;
}

// Differencing matrix:
// M = [ 1, 0, 0, ...
//      -1, 1, 0, ...
//      0, -1, 1, ...
//      ...      -1, 1 ]
// M is N x N
inline Eigen::MatrixXd differencingMatrix(int N)
{
  if (N == 0)
    throw std::invalid_argument("[linalg_utils::differencingMatrix] Size of specified matrix N is zero! N needs to be at least 1.");

  Eigen::MatrixXd M(N, N);
  M.setZero(N, N);

  for (int i = 0; i < N-1; i++)
  {
    M(i+1, i) = -1;
    M(i+1, i+1) = 1;
  }
  M(0, 0) = 1;
  M(N, N-1) = -1;

  return M;
}


// template <typename Derived,
//           template <typename, typename = std::allocator<Eigen::MatrixBase<Derived>>> class Container>
// inline Container<Eigen::MatrixBase<Derived>> Diff(const Container< Eigen::MatrixBase<Derived>>& vec)
// {
//   std::cout << "vec size: " << vec.size() << std::endl;
//   Container<Eigen::MatrixBase<Derived>> difference;
//   for (size_t i = 1; i < vec.size(); i++)
//   {
//     difference.push_back(vec[i]-vec[i-1]);
//   }
//
//   std::cout << "diff size: " << difference.size() << std::endl;
//   return difference;
// }

template <template <typename, typename = std::allocator<Eigen::Vector3d>> class Container>
inline Container<Eigen::Vector3d> Diff(const Container< Eigen::Vector3d>& vec)
{
  std::cout << "vec size: " << vec.size() << std::endl;
  Container<Eigen::Vector3d> difference;
  for (size_t i = 1; i < vec.size(); i++)
  {
    difference.push_back(vec[i]-vec[i-1]);
  }

  std::cout << "diff size: " << difference.size() << std::endl;
  return difference;
}

// template <typename Derived,
//           template <typename, typename = std::allocator<Eigen::MatrixBase<Derived>>> class Container>
// inline Container<Eigen::MatrixBase<Derived>> Dot(const Container< Eigen::MatrixBase<Derived>>& vec1, const Container< Eigen::MatrixBase<Derived>>& vec2)
// {
//   if (vec1.size() != vec2.size())
//   {
//     std::cout << "vec1 and vec2 are not the same size" << std::endl;
//   }
//
//   Container<Eigen::MatrixBase<Derived>> dot_product;
//   for (size_t i = 0; i < vec1.size(); i++)
//   {
//     dot_product.push_back(vec1[i].dot(vec2[i]));
//   }
//
//   return dot_product;
// }

// normalize each point one by one.
inline std::vector<Eigen::Vector3d> PointwiseNormalize(const std::vector<Eigen::Vector3d>& vec)
{
  std::vector<Eigen::Vector3d> vec_normalized;
  for (auto& p : vec)
  {
    vec_normalized.push_back(p.normalized());
  }
  return vec_normalized;
}

inline std::vector<double> Dot(const std::vector<Eigen::Vector3d>& vec1, const std::vector<Eigen::Vector3d>& vec2)
{
  if (vec1.size() != vec2.size())
  {
    std::cout << "vec1 and vec2 are not the same size" << std::endl;
  }

  std::vector<double> dot_product;
  for (size_t i = 0; i < vec1.size(); i++)
  {
    double dp = vec1[i].dot(vec2[i]);
    dot_product.push_back(dp);
  }

  return dot_product;
}

// Sum a vector
template<typename T>
T Sum(std::vector<T> &vec)
{
  T sum = 0;
  for (auto &elem : vec) {
    sum += elem;
  }
  return sum;
}


} // namespace linalg_utils
