#ifndef __CPP_UTILS_H__
#define __CPP_UTILS_H__

#include <Eigen/Core>
#include <iostream>
#include <random>
#include <vector>

namespace cpp_utils {
std::string methodName(const std::string &prettyFunction,
                       const std::string &file, const int &line);
std::string methodName(const std::string &prettyFunction);
template <typename T>
std::vector<T> arange(T start, T stop, T step = 1);
int select_random_element(const std::vector<int> &elements);
int select_random_element(const std::vector<int> &elements,
                          const std::vector<int> &weights);
int select_random_index(const std::vector<std::pair<int, int>> &elements);

int get_random_number(const int min, const int max);
std::vector<std::string> get_files_in_directory(const std::string &directory,
                                                const std::string &ext);
template <typename D>
void print_matrix(const Eigen::MatrixBase<D> &mat, const std::string &func,
                  const std::string &file, const int line);
}  // namespace cpp_utils

#endif  //__CPP_UTILS_H__
