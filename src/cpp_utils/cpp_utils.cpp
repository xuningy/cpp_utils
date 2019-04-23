#include "cpp_utils/cpp_utils.h"
#include <experimental/filesystem>

std::string cpp_utils::methodName(const std::string &prettyFunction)
{
  /**
    Purpose:  This function takes as input the "__PRETTY_FUNCTION__" macro
              from where its called and outputs the "class_name::function_name"

  */
  size_t begin, end;
  end = prettyFunction.find("(");
  begin = prettyFunction.substr(0, end).rfind(" ") + 1;
  end -= begin;
  return prettyFunction.substr(begin, end);
}

template std::vector<int>
cpp_utils::arange<int>(int start, int stop, int step);
template std::vector<float>
cpp_utils::arange<float>(float start, float stop, float step);
template std::vector<double>
cpp_utils::arange<double>(double start, double stop, double step);

template <typename T>
std::vector<T> cpp_utils::arange(T start, T stop, T step)
{
  std::vector<T> values;
  for (T value = start; value < stop; value += step)
    values.push_back(value);
  return values;
}

std::vector<std::string> cpp_utils::get_files_in_directory(
    const std::string &directory,
    const std::string &ext = "")
{
  std::vector<std::string> files;
  std::cout << "Size of Ext: " << ext.size() << std::endl;
  if (ext.size() == 0)
  {
    for (auto &p : std::experimental::filesystem::recursive_directory_iterator(directory))
    {
      files.push_back(p.path().string());
    }
  }
  else
  {
    for (auto &p : std::experimental::filesystem::recursive_directory_iterator(directory))
    {
      if (p.path().extension() == ext)
      {
        files.push_back(p.path().string());
      }
    }
  }

  return files;
}