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
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_set>

#include <Eigen/Core>
#include <Eigen/Geometry>

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
namespace stats_utils {

// Define a sufficiently small value, but
static constexpr float kSufficientlySmallFloat = 0.001;

template <typename T> using VecXt = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using MatXt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// RangeSample draws k samples of type T uniformly from the range [lb, ub), and
// returns an Eigen::VectorXd object. T is one of { int, float, double }.
template <typename T>
VecXt<T> RangeSample(const T lb, const T ub, const int k) {

  // Generate uniform real distribution from lb to ub.
  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_real_distribution<> distribution(lb, ub);

  // Sample from range.
  auto samples = VecXt<T>(k);
  for (int i = 0; i < k; i++) {
    samples[i] = distribution(generator);
  }
  return samples;
};

// RangeSample draws k samples of type T uniformly from the range [lb, ub), and
// returns a std::vector object. T is one of { int, float, double }.
template <typename T>
std::vector<T> RangeSample(const T lb, const T ub, const int k) {

  // Generate uniform real distribution from lb to ub.
  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_real_distribution<> distribution(lb, ub);

  // Sample from range.
  std::vector<T> samples;
  for (int i = 0; i < k; i++) {
    samples.push_back(distribution(generator));
  }
  return samples;
};

// DataSample draws k samples sampled uniformly at random from Container
// `data`, with replacement, and returns the sampled values in the same type
// Container.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
Container<T> DataSample(const Container<T> &data, const int k) {
  size_t N = data.size();

  // Check that the input arguments are valid.
  if (N == 0)
    throw std::invalid_argument("[stats_utils::DataSample] probability vector size needs to be larger than 1!");
  if (k <= 0)
    throw std::invalid_argument("[stats_utils::DataSample] Number to sample must be larger than 0!");

  // Generate uniform integer distribution.
  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_int_distribution<int> distribution(0, N - 1);

  // Sample from data.
  Container<T> sampled_data;
  for (int i = 0; i < k; i++) {
    int j = distribution(generator);
    sampled_data.push_back(data[j]);
  }
  return sampled_data;
}

// UniformDiscreteSample draws k samples sampled from a discrete range from 0
// to N, without replacement. Adapted from: https://stackoverflow.com/questions/28287138/c-randomly-sample-k-numbers-from-range-0n-1-n-k-without-replacement
// Originally from Robert Floyd  http://www.nowherenearithaca.com/2013/05/robert-floyds-tiny-and-beautiful.html
std::vector<int> UniformDiscreteSample(int N, int k) {
  std::random_device rd;
  std::mt19937 gen(rd());

  std::unordered_set<int> elems;
  for (int r = N - k; r < N; ++r) {
    int v = std::uniform_int_distribution<>(1, r)(gen);

    // there are two cases.
    // v is not in candidates ==> add it
    // v is in candidates ==> well, r is definitely not, because
    // this is the first iteration in the loop that we could've
    // picked something that big.

    if (!elems.insert(v).second) {
        elems.insert(r);
    }
  }

  // ok, now we have a set of k elements. but now
  // it's in a [unknown] deterministic order.
  // so we have to shuffle it:

  std::vector<int> result(elems.begin(), elems.end());
  std::shuffle(result.begin(), result.end(), gen);
  return result;
}

// DiscreteSample draws k samples sampled at random according to distribution
// `prob` with replacement, where `prob` is a probability array whose elements
// sum to 1. T is one of { float, double }.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
Container<int, std::allocator<int>> DiscreteSample(const Container<T> &prob,
                                                   const int k) {
  size_t N = prob.size();

  // Check that the input arguments are valid.
  if (N == 0)
    throw std::invalid_argument("[stats_utils::DiscreteSample] probability vector size needs to be larger than 1!");
  if (k <= 0)
    throw std::invalid_argument("[stats_utils::DiscreteSample] Number to sample must be larger than 0!");

  // Check that the sum of the probabilities is equal to 1 (accommodating
  // rounding errors).
  T prob_sum = 0;
  for (const T &p : prob)
    prob_sum += p;

  // TODO@Xuning: instead of throwing an error, normalize the vector instead.
  if (std::abs(prob_sum - 1.0) >= kSufficientlySmallFloat)
    throw std::invalid_argument("[stats_utils::DiscreteSample] Probability vector does not add up to 1!");

  // Create index vector.
  Container<int, std::allocator<int>> idx;
  for (size_t i = 0; i <= N; i += 1) {
    idx.push_back(i);
  }

  // Generate piecewise constant distribution.
  std::random_device rd;
  std::mt19937 generator(rd());
  std::piecewise_constant_distribution<> distribution(idx.begin(), idx.end(),
                                                      prob.begin());

  // Sample according to the probabilities.
  Container<int, std::allocator<int>> sampled_idx;
  for (int i = 0; i < k; i += 1) {
    sampled_idx.push_back(std::floor(distribution(generator)));
  }
  return sampled_idx;
}

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
  assert(samples.rows() == sample_dim);

  MatXt<T> centered = samples.colwise() - mean;
  MatXt<T> cov = (centered * centered.transpose()) / N;

  // Check that the sizes are correct.
  assert(cov.rows() == sample_dim);
  assert(cov.cols() == sample_dim);
  return cov;
}

// Covariance computes the weighted covariance for an Eigen Matrix `samples` of
// size MxN according to the mean `mean`, and outputs an Eigen matrix of size
// MxM. `samples` contains N samples of size M.
template <typename T>
MatXt<T> Covariance(const MatXt<T> &samples, const VecXt<T> &mean,
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

} // namespace stats_utils
