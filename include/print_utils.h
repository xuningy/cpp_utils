// Copyright 2018 Xuning Yang
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include <chrono>
#include <iostream>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace print_utils {

  inline void print(const std::vector<size_t> &vec, const std::string name = "") {
    std::cout << name << ": ";
    for (auto &v : vec) std::cout << v << " ";
    std::cout << std::endl;
  }

  inline void print(const std::vector<int> &vec, const std::string name = "") {
    std::cout << name << ": ";
    for (auto &v : vec) std::cout << v << " ";
    std::cout << std::endl;
  }

  inline void print(const std::vector<float> &vec, const std::string name = "") {
    std::cout << name << ": ";
    for (auto &v : vec) std::cout << v << " ";
    std::cout << std::endl;
  }

  inline void print(const std::vector<double> &vec, const std::string name = "") {
    std::cout << name << ": ";
    for (auto &v : vec) std::cout << v << " ";
    std::cout << std::endl;
  }

  template <typename Derived>
  inline void print(const Eigen::DenseBase<Derived> &mat, const std::string name = "") {
    Eigen::IOFormat OctaveFmt(4, 0, ", ", ";\n", "", "", "[", "]");
    Eigen::IOFormat OctaveVecFmt(4, 0, ", ", " ", "", "", "[", "]");

    std::cout << name << "(" << mat.rows() << "x" << mat.cols() << "): ";
    if (mat.rows() == 1 || mat.cols() == 1) std::cout << mat.format(OctaveVecFmt) << std::endl;
    else {
      std::cout << std::endl;
      std::cout << mat.format(OctaveFmt) << std::endl;
    }
  }

  inline void print(std::chrono::high_resolution_clock::time_point t1, std::chrono::high_resolution_clock::time_point t2) {
    std::chrono::duration<float> dur = t1 - t2;
    printf("(%.4fs)", dur.count());
  }

}
