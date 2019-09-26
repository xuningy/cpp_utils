/*
stats_utils.h
Copyright (C) 2019 Xuning Yang

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <cmath>
#include <chrono>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_set>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cpp_utils/linalg_utils.h>

// Defines a collection of statistics utility functions.
//
// Supported container types: Should be supported for all Sequence containers
// https://en.cppreference.com/w/cpp/container#Sequence_containers
//    std::array
//    std::vector
//    std::deque
//    std::forward_list
//    std::list
//
// Certain functions support only Eigen types only; see signatures below.
//

namespace lu = linalg_utils;

namespace stats_utils {

// Define a sufficiently small value, but
static constexpr float kSufficientlySmallFloat = 0.01;

template <typename T> using VecXt = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using MatXt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// GaussianPdf computes the probability density for a sample `x`, according to
// mean `mean` and covariance `sigma`.
// Inputs:
//  x:      vector of size M
//  mean:   vector of size M
//  sigma:  matrix of size MxM
template <typename T>
T GaussianPdf(const VecXt<T> &x, const VecXt<T> &mean, const MatXt<T> &sigma) {
  if (x.rows() != mean.rows() || x.cols() != mean.cols())
    throw std::invalid_argument("[stats_utils::GaussianPdf] Size of x does not equal to size of mean (should be a vector of n elements each)");
  if (sigma.rows() != sigma.cols())
    throw std::invalid_argument("[stats_utils::GaussianPdf] Sigma is not square.");

  // Precompute constants.
  int n = std::max(x.rows(), x.cols());
  const T log_sqrt_2pi = 0.5 * std::log(2 * M_PI);

  // Use Cholesky decomposition so as to prevent numerical instability in
  // computing the inverse and determinant.
  Eigen::LLT<MatXt<T>> chol(sigma);

  if (chol.info() == Eigen::NumericalIssue)
    throw std::runtime_error("Possibly non PSD matrix!");

  const typename Eigen::LLT<MatXt<T>>::Traits::MatrixL &L = chol.matrixL();
  T dist = (L.solve(x - mean)).squaredNorm();
  return std::exp(-n * log_sqrt_2pi - 0.5 * dist) / L.determinant();
}

// Average computes the average and stddev of the vector.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
T Average(const Container<T>& v, T *stddev) {
  if (v.empty())
  {
    std::cout << "[stats_utils::Average] WARNING: container has no elements; returning NaN." << std::endl;
    *stddev = std::numeric_limits<T>::quiet_NaN();
    return std::numeric_limits<T>::quiet_NaN();
  }

  T sum = std::accumulate(v.begin(), v.end(), 0.0);
  T mean = sum / v.size();

  Container<T> diff(v.size());
  std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
  T sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  *stddev = std::sqrt(sq_sum / v.size());

  return mean;
}

// WeightedMean computes the weighted mean of samples of size MxN according to
// weights of size N, and outputs a vector mean of size M.
template <typename T> VecXt<T> WeightedMean(const MatXt<T> &samples, const VecXt<T> &weights)
{
  int N = samples.cols();
  if (N != weights.size())
    throw std::invalid_argument(std::string("[stats_utils::WeightedMean] number of samples and number of weights are not the same: %d vs. %d!", N, weights.size()));

  T weights_sum = weights.sum();

  // Compute the weighted samples.
  MatXt<T> weighted_samples = samples.array().rowwise() * weights.transpose().array();

  // Compute the mean.
  VecXt<T> weighted_mean = 1.0 / weights_sum * weighted_samples.rowwise().sum();

  return weighted_mean;
}

// Covariance computes the sample covariance for an Eigen Matrix `samples` of
// size MxN, and outputs an Eigen matrix of size MxM.
// `samples` contains N samples of size M.
template <typename T> MatXt<T> Covariance(const MatXt<T> &samples) {
  int N = samples.cols();
  VecXt<T> mean = samples.rowwise().mean();
  int sample_dim = mean.size();
  MatXt<T> centered = samples.colwise() - mean;
  MatXt<T> cov = (centered * centered.transpose()) / N;

  return cov;
}

// Covariance computes the sample covariance for an Eigen Matrix `samples` of
// size MxN according to the mean `mean`, and outputs an Eigen matrix of size
// MxM. `samples` contains N samples of size M.
template <typename T>
MatXt<T> Covariance(const MatXt<T> &samples, const VecXt<T> &mean) {
  int N = samples.cols();
  int sample_dim = mean.size();
  if (samples.rows() != sample_dim) {
    throw std::invalid_argument("[stats_utils::Covariance] sample data has size %d while mean has data size %d. Exiting.", samples.rows(), sample_dim);
  }

  MatXt<T> centered = samples.colwise() - mean;
  MatXt<T> cov = (centered * centered.transpose()) / N;

  // Check that the sizes are correct.
  assert(cov.rows() == sample_dim);
  assert(cov.cols() == sample_dim);
  return cov;
}

// WeightedCovariance computes the weighted covariance for an Eigen Matrix `samples` of
// size MxN according to the mean `mean`, and outputs an Eigen matrix of size
// MxM. `samples` contains N samples of size M.
template <typename T>
MatXt<T> WeightedCovariance(const MatXt<T> &samples, const VecXt<T> &mean,
                    const VecXt<T> &weights) {

  int N = samples.cols();
  int sample_dim = mean.size();
  assert(samples.rows() == sample_dim);
  assert(weights.size() == N);

  MatXt<T> centered = samples.colwise() - mean;
  MatXt<T> weighted_centered = centered.array().rowwise() * weights.transpose().array();
  MatXt<T> cov = (centered * weighted_centered.transpose()) / weights.sum();

  // Check that the sizes are correct.
  assert(cov.rows() == sample_dim);
  assert(cov.cols() == sample_dim);
  return cov;
}

// PercentileValue computes the value at the specified percentile.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
T PercentileValue(const Container<T> &v, float percentile) {

  if (percentile < 0 || percentile > 1.0) {
    throw std::invalid_argument("[stats_utils::PercentileValue] percentile needs to be between 0.0 and 1.0!");
  }
  Container<T, std::allocator<T>> v_sorted = v;
  std::sort(v_sorted.begin(), v_sorted.end());

  if (v.size() <= 2) {
    std::cout << "[WARNING][stats_utils::PercentileValue] Vector too short. Returning the largest value." << std::endl;
    return v_sorted.back();
  }

  // compute the index value
  T R = percentile * (v_sorted.size() + 1);

  // interpolate
  unsigned int IR = std::max((int)std::floor(R), 1);
  T FR = R - (int)std::floor(R);

  size_t LB = IR - 1;
  size_t UB = IR;

  return v_sorted[LB] + (v_sorted[UB] - v_sorted[LB]) * FR;

}

// KolmogorovSmirnov performs the Kolmogorov-Smirnov test between two discrete distributions f & g. KS first computes the CDF for each, and then computes the maximum distance.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
T KolmogorovSmirnov(const Container<T>& f, const Container<T>& g) {
  if (f.size() != g.size()) {
    throw std::invalid_argument(std::string("[stats_utils::KolmogorovSmirnov] the discrete distributions f and g have different sizes: f.size() = %d, g.size() = %d", f.size(), g.size()));
  }

  size_t N = f.size();
  T dx = 1 / (T)N;
  Container<T> F = lu::Cumtrapz(dx, f);
  Container<T> G = lu::Cumtrapz(dx, g);
  T sumf = 0.0;
  for (auto & elem : F) {
    sumf += elem;
  }

  T sumg = 0.0;
  for (auto & elem : G) {
    sumg += elem;
  }

  if ( std::abs(sumf - 1.0) > 0.01 || std::abs(sumg - 1.0) > 0.01) {
    throw std::invalid_argument("[stats_utils::KolmogorovSmirnov] f or g is an improper distribution!");
  }

  // Find the max vertical distance between the two values.
  T max_dist = 0;
  for (size_t i = 0; i < N; i++) {
    T dist = std::abs(F[i] - G[i]);

    std::cout << F[i] << " - " << G[i] << " = " << dist << std::endl;
    if (dist > max_dist) {
      max_dist = dist;
    }
  }

  return max_dist;
}

} // namespace stats_utils
