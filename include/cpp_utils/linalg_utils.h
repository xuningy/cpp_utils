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
  }

  // Check that the input arguments are valid.
  if (N == 0)
    throw std::invalid_argument("[linalg_utils::Linspace] probability vector size needs to be larger than 1!");
  if (lb >= ub)
    throw std::invalid_argument("[linalg_utils::Linspace] Number to sample must be larger than 0!");

  // Generate the evenly spaced vector.
  T h = (ub - lb) / static_cast<T>(N - 1);
  std::vector<T> vec;
  for (int n = 0; n < N; n++) {
    vec.push_back(lb + n*h);
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


} // namespace linalg_utils
