/*
sample_utils.h
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

// Defines a collection of sample utility functions.
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
namespace sample_utils {

// Define a sufficiently small value, but
static constexpr float kSufficientlySmallFloat = 0.01;

template <typename T> using VecXt = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using MatXt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// RangeSample draws k samples of type T uniformly from the range [lb, ub), and
// returns an Eigen::VectorXd object. T is one of { int, float, double }.

template <typename T>
VecXt<T> RangeSampleEigen(const T lb, const T ub, const int k, std::mt19937& generator) {

  // Generate uniform real distribution from lb to ub.
  std::uniform_real_distribution<> distribution(lb, ub);

  // Sample from range.
  auto samples = VecXt<T>(k);
  for (int i = 0; i < k; i++) {
    samples[i] = distribution(generator);
  }
  return samples;
};

template <typename T>
VecXt<T> RangeSampleEigen(const T lb, const T ub, const int k) {

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);
  return RangeSample(lb, ub, k, generator);
};

// RangeSample draws k samples of type T uniformly from the range [lb, ub), and
// returns a std::vector object. T is one of { int, float, double }.

template <typename T>
std::vector<T> RangeSample(const T lb, const T ub, const int k, std::mt19937& generator) {

  // Generate uniform real distribution from lb to ub.
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
Container<T> DataSample(const Container<T> &data, const int k, std::mt19937& generator) {
  size_t N = data.size();

  // Check that the input arguments are valid.
  if (N == 0)
    throw std::invalid_argument("[sample_utils::DataSample] Probability vector is empty!");
  if (k <= 0)
    throw std::invalid_argument("[sample_utils::DataSample] Number to sample < 0!");

  // Generate uniform integer distribution.
  std::uniform_int_distribution<int> distribution(0, N - 1);

  // Sample from data.
  Container<T> sampled_data;
  for (int i = 0; i < k; i++) {
    int j = distribution(generator);
    sampled_data.push_back(data[j]);
  }
  return sampled_data;
}

template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
Container<T> DataSample(const Container<T> &data, const int k) {
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);
  return DataSample(data, k, generator);
}

// UniformDiscreteSample draws k samples sampled from a discrete range from 0
// to N, without replacement. Adapted from: https://stackoverflow.com/questions/28287138/c-randomly-sample-k-numbers-from-range-0n-1-n-k-without-replacement
// Originally from Robert Floyd  http://www.nowherenearithaca.com/2013/05/robert-floyds-tiny-and-beautiful.html
inline std::vector<int> UniformDiscreteSample(int N, int k, std::mt19937& generator) {

  std::unordered_set<int> elems;
  for (int r = N - k; r < N; ++r) {
    int v = std::uniform_int_distribution<>(1, r)(generator);

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
  std::shuffle(result.begin(), result.end(), generator);
  return result;
}


inline std::vector<int> UniformDiscreteSample(int N, int k) {
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);
  return UniformDiscreteSample(N, k, generator);
}

// DiscreteSample draws k samples sampled at random according to
// distribution `prob` with replacement, where `prob` is a probability array
// whose elements sum to 1. T is one of { float, double }.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
Container<int, std::allocator<int>> DiscreteSample(
  const Container<T> &prob, const int k, std::mt19937& generator) {
  size_t N = prob.size();

  // Check that the input arguments are valid.
  if (N == 0)
    throw std::invalid_argument("[sample_utils::DiscreteSample] probability vector size needs to be larger than 1!");
  if (k <= 0)
    throw std::invalid_argument("[sample_utils::DiscreteSample] Number to sample must be larger than 0!");

  // Check that the sum of the probabilities is equal to 1 (accommodating
  // rounding errors).
  T prob_sum = 0;
  for (const T &p : prob)
    prob_sum += p;

  // TODO@Xuning: instead of throwing an error, normalize the vector instead.
  if (std::abs(prob_sum - 1.0) >= kSufficientlySmallFloat) {
    std::cout << "Probability vector sum: " << prob_sum << std::endl;
    throw std::invalid_argument("[sample_utils::DiscreteSample] Probability vector does not add up to 1!");
  }


  // Create index vector.
  Container<int, std::allocator<int>> idx;
  for (size_t i = 0; i < N; i += 1) {
    idx.push_back(i);
  }

  // Generate piecewise constant distribution.
  std::piecewise_constant_distribution<> distribution(idx.begin(), idx.end(),
                                                      prob.begin());

  // Sample according to the probabilities.
  Container<int, std::allocator<int>> sampled_idx;
  for (int i = 0; i < k; i += 1) {
    sampled_idx.push_back(std::floor(distribution(generator)));
  }
  return sampled_idx;
}

template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
Container<int, std::allocator<int>> DiscreteSample(
  const Container<T> &prob, const int k) {
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);

    return DiscreteSample(prob, k, generator);
  }


// DiscreteSampleWithoutReplacement draws k samples sampled at random according to distribution
// `prob` without replacement, where `prob` is a probability array whose
// elements sum to 1. T is one of { float, double }.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
Container<int, std::allocator<int>> DiscreteSampleWithoutReplacement(
  const Container<T> &prob, const int k, std::mt19937& generator) {
  size_t N = prob.size();

  // Check that the input arguments are valid.
  if (N == 0)
    throw std::invalid_argument("[sample_utils::DiscreteSampleWithoutReplacement] probability vector size needs to be larger than 1!");
  if (k <= 0)
    throw std::invalid_argument("[sample_utils::DiscreteSampleWithoutReplacement] Number to sample must be larger than 0!");

  // Check that the sum of the probabilities is equal to 1 (accommodating
  // rounding errors).
  T prob_sum = 0;
  for (const T &p : prob)
    prob_sum += p;

  // TODO@Xuning: instead of throwing an error, normalize the vector instead.
  if (std::abs(prob_sum - 1.0) >= kSufficientlySmallFloat) {
    std::cout << "Probability vector sum: " << prob_sum << std::endl;
    throw std::invalid_argument("[sample_utils::DiscreteSampleWithoutReplacement] Probability vector does not add up to 1!");
  }


  // Create index vector.
  Container<int, std::allocator<int>> idx;
  for (size_t i = 0; i < N; i++) {
    idx.push_back(i);
  }

  // Generate discrete distribution.
  Container<int, std::allocator<int>> sampled_idx;

  // If the number of sampled elements is greater than the values themselves,
  // just shuffle (uniformly random, not weighted) the index vector and return.
  if (k >= N) {

    sampled_idx = idx;
    std::shuffle(sampled_idx.begin(), sampled_idx.end(), generator);
    return sampled_idx;
  }

  // Sample according to the probabilities.
  std::discrete_distribution<> distribution(prob.begin(), prob.end());

  while (sampled_idx.size() < k) {
    int value = distribution(generator);
    if (std::find(sampled_idx.begin(), sampled_idx.end(), value) == sampled_idx.end()) {
      sampled_idx.push_back(value);
    }
  }
  return sampled_idx;
}


template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
Container<int, std::allocator<int>> DiscreteSampleWithoutReplacement(
  const Container<T> &prob, const int k) {
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);

    return DiscreteSampleWithoutReplacement(prob, k, generator);
  }


  // std::vector<unsigned int> UserModel::sampleWithReplacement(const unsigned int K, const unsigned int N)
  // {
  //   std::vector<unsigned int> sampled_idx;
  //
  //   // importance sampling with replacement (includes duplicates in sampled_idx)
  //
  //   std::srand(std::time(0));
  //   for (unsigned int j = 0; j < K; j++) {
  //     double total = wgt.sum();
  //     double val = (double)std::rand()/(double)RAND_MAX*total;
  //     for (unsigned int k= 0; k < N; k++)
  //     {
  //       val = val - wgt(k);
  //       if (val <= 0) {
  //         sampled_idx.push_back(k);
  //         break;
  //       }
  //     }
  //
  //   }
  //
  //   // add 5 fixed motion primitives for N = 101
  //   // sampled_idx.push_back(0);
  //   // sampled_idx.push_back(101);
  //   // sampled_idx.push_back(51);
  //   // sampled_idx.push_back(26);
  //   // sampled_idx.push_back(77);
  //
  //   return sampled_idx;
  //
  // }
  //
  // std::vector<unsigned int> UserModel::sampleWithoutReplacement(const unsigned int K, const unsigned int N)
  // {
  //   std::vector<unsigned int> sampled_idx;
  //
  //   Eigen::VectorXd wgt2 = wgt;
  //
  //   // importance sampling without replacement (no duplicates in sampled_idx)
  //
  //   std::srand(std::time(0));
  //   for (unsigned int j = 0; j < K; j++) {
  //     double total = wgt2.sum();
  //     double val = (double)std::rand()/(double)RAND_MAX*total;
  //     for (unsigned int k= 0; k < N; k++)
  //     {
  //       val = val - wgt2(k);
  //       if (val <= 0) {
  //         wgt2(k) = 0;
  //         sampled_idx.push_back(k);
  //         break;
  //       }
  //     }
  //
  //   }
  //
  //   // sanity check
  //
  //   if (sampled_idx.size() != K) ROS_ERROR_STREAM("[sampleWithoutReplacement] length of sampled_idx " << sampled_idx.size() << " isn't the same as K " << K);
  //
  //   return sampled_idx;
  //
  // }
  //
  // std::vector<unsigned int> UserModel::sampleFromInverseCDF(const unsigned int K, const unsigned int N)
  // {
  //   std::vector<unsigned int> sampled_idx;
  //
  //   Eigen::VectorXd x;
  //   x.setLinSpaced(N, 0, 1);
  //
  //   Eigen::VectorXd f_unif;
  //   f_unif.setLinSpaced(K, 0, 1);
  //
  //   double dx = x(1) - x(0);
  //   Eigen::VectorXd cdf = lu::Cumtrapz(dx, wgt);
  //
  //   cdf = cdf.array()/cdf.maxCoeff();
  //
  //   for (unsigned int k = 0; k < K; k++)
  //   {
  //     Eigen::VectorXd unif_cdf(N);
  //     unif_cdf.setConstant(f_unif(k));
  //     unsigned int idx;
  //     (cdf - unif_cdf).array().abs().matrix().minCoeff(&idx);
  //
  //     sampled_idx.push_back(idx);
  //   }
  //
  //   return sampled_idx;
  //
  // }

} // namespace sample_utils
