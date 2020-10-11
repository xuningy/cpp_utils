/*
vector_utils.h
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

#include <algorithm>
#include <deque>
#include <numeric>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

template <typename T> using VecXt = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using MatXt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// Defines a collection of vector utility functions.
//
// Supported container types: Should be supported for all Sequence containers
// https://en.cppreference.com/w/cpp/container#Sequence_containers
//    std::array
//    std::vector
//    std::deque
//    std::forward_list
//    std::list
//

namespace vector_utils {

// Sorts a VecContainer `v` of type `T` in ascending order, and outputs the
// sorted result `v_sorted` in the same container type VecContainer, and the
// sorted indices `idx_sorted` in container type IdxContainer.
template <typename T,
          template <typename, typename = std::allocator<T>> class VecContainer,
          template <typename, typename = std::allocator<size_t>> class IdxContainer>
void Sort(const VecContainer<T> &v, VecContainer<T> *v_sorted,
          IdxContainer<size_t> *idx_sorted) {

  // Initialize original index locations.
  IdxContainer<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  // Sort indexes based on comparing values in v.
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  // Populate the output parameters.
  *idx_sorted = idx;
  *v_sorted = v;
  std::sort(v_sorted->begin(), v_sorted->end());
}

template <typename T,
          template <typename, typename = std::allocator<T>> class VecContainer,
          template <typename, typename = std::allocator<float>> class ArgContainer>
void ArgSort(const VecContainer<T> &v, VecContainer<T> *v_sorted, const ArgContainer<float> &varg, ArgContainer<float> *varg_sorted) {
  std::vector<size_t> args_sorted_idx;
  Sort<float>(varg, varg_sorted, &args_sorted_idx);
  // Push escape points by sorted cost
  for(uint i = 0; i < v.size(); i++) {
    v_sorted->push_back(v[args_sorted_idx[i]]);
  }
}

// =========================== FIND OPERATIONS =============================

// Find an element in the Container, and return it's position
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
bool Find(const Container<T> &v, const T &elem, size_t *index) {

  for (auto iter = v.begin(); iter != v.end(); ++iter) {
    if (*iter == elem) {
      *index = iter - v.begin();
      return true;
    }
  }
  return false;
}

// Find an element in the Container.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
bool Find(const Container<T> &v, const T &elem) {
  size_t index;
  return Find(v, elem, &index);
}

// Find the indices of all the elements that are equivalent to elem.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
Container<int> FindIndices(const Container<T> &v, const T &elem) {
  Container<int> locations;
    for (auto iter = v.begin(); iter != v.end(); ++iter) {
      if (*iter == elem) {
        locations.push_back(iter - v.begin());
      }
    }
  return locations;
}

// Find the indices of all the elements that satisfies the predicate.
template <typename T,
          class UnaryPredicate,
          template <typename, typename = std::allocator<T>> class Container>
Container<int> Find(const Container<T> &v, UnaryPredicate* predicate) {
  Container<int> locations;
    for (auto iter = v.begin(); iter != v.end(); ++iter) {
      if (predicate(*iter)) {
        locations.push_back(iter - v.begin());
      }
    }
  return locations;
}

// =========================== VALUE OPERATIONS  =============================

// Find max element in a Container v, and return it's value and location.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
T Max(const Container<T> &v, size_t *location) {
  Container<T, std::allocator<T>> v_sorted;
  Container<size_t, std::allocator<size_t>> idx_sorted;

  Sort(v, &v_sorted, &idx_sorted);
  *location = idx_sorted.back();
  return v_sorted.back();
}

// Find max element in a Container v, and return it's value.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
T Max(const Container<T> &v) {
  return *std::max_element(v.begin(), v.end());
}

// Find min element in a Container v, and return it's value and location.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
T Min(const Container<T> &v, size_t *location) {
  Container<T, std::allocator<T>> v_sorted;
  Container<size_t, std::allocator<size_t>> idx_sorted;

  Sort(v, &v_sorted, &idx_sorted);
  *location = idx_sorted.front();
  return v_sorted.front();
}

// Find min element in a Container v, and return it's value.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
T Min(const Container<T> &v) {
  return *std::min_element(v.begin(), v.end());
}

// Find min and max element in a Container v, and return it's values as a
// std::tuple<T, T>{min, max}.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
std::tuple<T, T> MinMax(const Container<T> &v, std::tuple<T, T> *locations) {
  Container<T, std::allocator<T>> v_sorted;
  Container<size_t, std::allocator<size_t>> idx_sorted;

  Sort(v, &v_sorted, &idx_sorted);
  *locations = std::make_tuple(idx_sorted.front(), idx_sorted.back());
  return std::make_tuple(v_sorted.front(), v_sorted.back());
}

// Find min and max element in a Container v, and return it's values as a
// std::tuple<T, T>{min, max}.
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
std::tuple<T, T> MinMax(const Container<T> &v) {
  std::tuple<T, T> locations;
  return MinMax(v, &locations);
}

// Invert all the values in `vec` with type InContainer, wrt to the max value
// in `vec`; i.e. (max(vec) - vec). The inverted output is stored in
// `inverted_vec`, with container type OutContainer.
template <typename T,
          template <typename, typename = std::allocator<T>> class InContainer,
          template <typename, typename = std::allocator<T>> class OutContainer>
void Invert(const InContainer<T>& vec, OutContainer<T> *inverted_vec) {
  T max = Max(vec);

  inverted_vec->clear();
  for (const T &elem : vec) {
    inverted_vec->push_back(max - elem);
  }

  return;
}

// Normalize all the values in `vec` with type InContainer, wrt to the total
// value in `vec`; i.e. (vec/sum). The output is stored in `normalized_vec`,
// with container type OutContainer.
template <typename T,
          template <typename, typename = std::allocator<T>> class InContainer,
          template <typename, typename = std::allocator<T>> class OutContainer>
void Normalize(const InContainer<T>& vec, OutContainer<T> *normalized_vec) {
  T sum = 0;
  for (auto &elem : vec) {
    sum += elem;
  }
  normalized_vec->clear();
  if (sum == 1) {
    *normalized_vec = vec;
  } else if (sum == 0) {
    normalized_vec->assign(vec.size(), 1.0/vec.size());
  } else {
    T factor = 1.0/sum; // doing a single division and then multiplication lowers the numerical instability associated with division.
    for (const T &elem : vec) {
      normalized_vec->push_back(elem * factor);
    }
  }

  return;
}

// RescaleMinMax rescales all the values in `vec` with type InContainer wrt min
// max value in `vec` such that all values are bounded between 0 and 1. The
// output is stored in `rescaled_vec`, with container type OutContainer.
template <typename T,
          template <typename, typename = std::allocator<T>> class InContainer,
          template <typename, typename = std::allocator<T>> class OutContainer>
void RescaleMinMax(const InContainer<T>& vec, OutContainer<T> *rescaled_vec) {
  std::tuple<T, T> minmax = MinMax(vec);
  T min = std::get<0>(minmax);
  T max = std::get<1>(minmax);

  rescaled_vec->clear();
  if ((min == 0 && max == 1) || min == max) {
    *rescaled_vec = vec;
  } else {
    T factor = 1.0/(max - min); // doing a single division and then multiplication lowers the numerical instability associated with division.
    for (const T &elem : vec) {
      rescaled_vec->push_back((elem - min) * factor);
    }
  }

  return;
}

// RescaleMinMax rescales all the values in `v`
template <typename T>
VecXt<T> RescaleMinMax(const VecXt<T>& v) {
  T min = v.minCoeff();
  T max = v.maxCoeff();

  if ((min == 0 && max == 1) || max == min) return v;

  return (v.array() - VecXt<T>::Constant(v.size(), 1, min).array() ) / (max - min);
}

// =========================== REMOVE OPERATIONS =============================

// Remove all instances of a value from a vector
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
void RemoveAll(Container<T>& vec, T val)
{
  vec.erase(std::remove(vec.begin(), vec.end(), val), vec.end());
}

// Remove all elements that meet the unary predicate's condition
template <typename T,
          class UnaryPredicate,
          template <typename, typename = std::allocator<T>> class Container>
void RemoveAllIf(Container<T>& vec, UnaryPredicate predicate)
{
  vec.erase(std::remove_if(vec.begin(), vec.end(), predicate), vec.end());
}

// Remove one index from a vector
template <typename T,
          template <typename, typename = std::allocator<T>> class Container>
void RemoveAtIndex(Container<T>& vec, int index)
{
  vec.erase(vec.begin() + index);
}

// Remove a list of indices from a vector by copying data to a new vector in blocks.
template <typename T,
          template <typename, typename = std::allocator<int>> class Container>
std::vector<T> RemoveAtIndices(const std::vector<T>& v1, Container<int>& indices)
{
  // if indices empty, do nothing.
  if (indices.empty()) return v1;

  // Create a new vector
  std::vector<T> v2;
  v2.reserve(v1.size() - indices.size());

  // Sort indices.
  std::sort(indices.begin(), indices.end());

  // Copy blocks over at once.
  typename std::vector<T>::const_iterator itBlockBegin = v1.begin();
  for (typename Container<int>::const_iterator it = indices.begin(); it != indices.end(); ++it)
  {
    typename std::vector<T>::const_iterator itBlockEnd = v1.begin() + *it;
    if (itBlockBegin != itBlockEnd)
    {
      std::copy(itBlockBegin, itBlockEnd, std::back_inserter(v2));
    }
    itBlockBegin = itBlockEnd + 1;
  }

  // Copy last block.
  if(itBlockBegin != v1.end())
  {
    std::copy(itBlockBegin, v1.end(), std::back_inserter(v2));
  }

  return v2;
}

// Remove a list of indices from a deque by iterating over the deque.
template <typename T,
          template <typename, typename = std::allocator<int>> class Container>
void RemoveAtIndices(std::deque<T>& v, Container<int>& indices)
{
  // If indices empty, do nothing.
  if (indices.empty()) return;

  // If only one index, directly remove it.
  if (indices.size() == 1)
  {
    RemoveAtIndex(v, indices[0]);
    return;
  }

  // Sort indices.
  std::sort(indices.begin(), indices.end());

  int i = 0;
  int index = indices[i];

  //  Iterate over deque and remove the appropriate elements.
  for (typename std::deque<T>::const_iterator it = v.begin(); it != v.end(); ) {
    if ( (it - v.begin() - i) == index) {
      it = v.erase(it); // returns the next iterator

      // update the next index to find.
      i++;
      if (i >= indices.size()) break;
      index = indices[i];
    } else {
      ++it; // increment the iterator
    }
  }

  return;
}

// Return a slice from position m to n inclusive. Credit: https://www.techiedelight.com/get-slice-sub-vector-from-vector-cpp/
template<typename T>
std::vector<T> Slice(std::vector<T> &v, int m, int n)
{
  std::vector<T>vec(n-m+1);
  std::copy(v.begin() + m, v.begin() + n + 1, vec.begin());
  return vec;
}

} // namespace vector_utils
