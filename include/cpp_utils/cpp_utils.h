#ifndef __CPP_UTILS_H__
#define __CPP_UTILS_H__

#include <iostream>
#include <vector>

namespace cpp_utils
{
std::string methodName(const std::string &prettyFunction);
template <typename T>
std::vector<T> arange(T start, T stop, T step = 1);
std::vector<std::string> get_files_in_directory(const std::string &directory,
                                                const std::string &ext);
} // namespace cpp_utils

#endif //__CPP_UTILS_H__