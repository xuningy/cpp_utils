/*
divergences.h
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

template <typename T> using VecXt = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using MatXt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

namespace divergence {

// KullbackLeiblerGaussian computes the Kullback Leibler divergence between two
// univariate Gaussian distributions.
template <typename T>
T KullbackLeiblerGaussian(T mean_a, T mean_b, T sigma_a, T sigma_b)
{
  T kld = std::log(sigma_b / sigma_a)
        + (sigma_a * sigma_a + (mean_a - mean_b) * (mean_a - mean_b)) / (2 * sigma_b * sigma_b) - 0.5;

  return kld;
}

// KullbackLeiblerGaussian computes the Kullback Leibler divergence between two
// multivariate Gaussian distributions.
template <typename T>
T KullbackLeiblerGaussian(const VecXt<T> &mean_a, const VecXt<T> &mean_b,
                  const MatXt<T> &cov_a, const MatXt<T> &cov_b)
{
  int d = mean_a.size();

  MatXt<T> cov_b_inv = cov_b.inverse();
  T kld = 0.5 * (std::log(cov_b.determinant() / cov_a.determinant())
                - d
                + (cov_b_inv * cov_a).trace()
                + ((mean_b - mean_a).transpose() * cov_b_inv * (mean_b - mean_a)));
  return kld;
}

// KullbackLeiblerGaussian computes the Kullback Leibler divergence given two
// pdf's.
template <typename T>
T KullbackLeibler(const VecXt<T> &p, const VecXt<T> &q)
{
  if (p.size() != q.size()) {
    throw std::invalid_argument(std::string("[divergence::KL] p and q pdf's are of different sizes: p.size() = %d, q.size() = %d", p.size(), q.size()));
  }
  // Normalize the two distributions
  VecXt<T> pn = p.normalized();
  VecXt<T> qn = q.normalized();

  return pn.array() * ((pn.array() / qn.array()).matrix().log()).array();
}

// JeffreysGaussian computes the Jeffreys divergence between two univariate
// Gaussian distributions.
template <typename T>
T JeffreysGaussian(T mean_a, T mean_b, T sigma_a, T sigma_b)
{
  return KullbackLeiblerGaussian(mean_a, mean_b, sigma_a, sigma_b)
         + KullbackLeiblerGaussian(mean_b, mean_a, sigma_b, sigma_a);
}

// JeffreysGaussian computes the Jeffreys divergence between two multivariate
// Gaussian distributions.
template <typename T>
T JeffreysGaussian(const VecXt<T> &mean_a, const VecXt<T> &mean_b,
                  const MatXt<T> &cov_a, const MatXt<T> &cov_b)
{
  return KullbackLeiblerGaussian(mean_a, mean_b, cov_a, cov_b)
         + KullbackLeiblerGaussian(mean_b, mean_a, cov_b, cov_a);
}

// JensenShannon computes the Jensen Shannon divergence between two pdf's.
template <typename T>
T JensenShannon(const VecXt<T> &p, const VecXt<T> &q)
{
  VecXt<T> m = (p + q) * 0.5;
  return 0.5 * (KullbackLeibler(p, m)
                + KullbackLeibler(q, m));
}

} // namespace divergence

} // namespace stats_utils
