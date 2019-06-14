#include "cpp_utils/cpp_utils.h"
#include <experimental/filesystem>
#include <random>
#define getName(var) #var

std::random_device rd_;
std::mt19937 gen_(rd_());

std::string cpp_utils::methodName(const std::string &prettyFunction,
                                  const std::string &file,
                                  const int &line)
{
  /**
    Purpose:  This function takes as input the "__PRETTY_FUNCTION__" macro
              from where its called and outputs the "class_name::function_name"

  */
  size_t begin, end;
  end = prettyFunction.find("(");
  begin = prettyFunction.substr(0, end).rfind(" ") + 1;
  end -= begin;

  size_t f_begin, f_end;
  f_end = file.size();
  f_begin = file.substr(0, file.size()).find_last_of("/") + 1;
  std::string ret = "File: ";
  ret.append(file.substr(f_begin, f_end));
  ret.append("\t Function: ");
  ret.append(prettyFunction.substr(begin, end));
  ret.append("\t Line No: ");
  ret.append(std::to_string(line));
  // std::string ret = prettyFunction.substr(begin, end);

  return ret;
}

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

template void cpp_utils::print_matrix<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(const Eigen::MatrixBase<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &mat, const std::string &func, const std::string &file, const int line);
template void cpp_utils::print_matrix<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(const Eigen::MatrixBase<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &mat, const std::string &func, const std::string &file, const int line);
template void cpp_utils::print_matrix<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(const Eigen::MatrixBase<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &mat, const std::string &func, const std::string &file, const int line);

template <typename D>
void cpp_utils::print_matrix(const Eigen::MatrixBase<D> &mat, const std::string &func,
                              const std::string &file, const int line)
{
  /**
  * @brief:  Print out a matrix, the function, the file and the line number it's called from
  * 
  * This method takes as input an Eigen matrix and prints it out in a numpy format along with 
  * the debug message of the file, function and line number the print out is called at.
  * 
  * @param  mat is the input matrix to be printed
  * @param func is the C++ macro __PRETTY_FUNCTION__
  * @param file is the C++ macro __FILE__
  * @param line is the C++ macro __LINE__
  * 
  * @return void
  */
  std::string print_message = cpp_utils::methodName(func, file, line);
  printf("[%s]: Matrix: %s \n", print_message.c_str(), getName(mat));
  std::cout << "[";
  for (int i = 0; i < mat.rows()-1; ++i)
  {
    if(i != 0)
      std::cout << " [ ";
    else 
      std::cout << "[ ";
    for (int j = 0; j < mat.cols() - 1; ++j)
    {
      std::cout << mat(i, j) << ", ";
    }
    std::cout << mat(i, mat.cols() - 1);
    std::cout << "]\n";
  }
  std::cout << " [ ";
  for (int j = 0; j < mat.cols() - 1; ++j)
  {
    std::cout << mat(mat.rows() - 1, j) << ", ";
  }
  std::cout << mat(mat.rows() - 1, mat.cols() - 1);
  std::cout << "]";
  std::cout << "]\n\n";
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

int cpp_utils::select_random_element(const std::vector<int> &el)
{
  std::uniform_int_distribution<> dis(0, (int)el.size() - 1);
  int index = dis(gen_);
  return index;
}

int cpp_utils::get_random_number(const int min, const int max)
{
  std::uniform_int_distribution<> dis(min, max - 1);
  int index = dis(gen_);
  return index;
}