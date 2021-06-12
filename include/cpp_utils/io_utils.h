/*
io_utils.h
Copyright (C) 2020 Xuning Yang

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

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace io_utils {

// gets a string for datetime in YYYY-mm-dd_hh-mm (24h) format. Check documentation for details: http://www.cplusplus.com/reference/ctime/strftime/
inline std::string DateTime()
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  // date in YYYY-mm-dd_HH-mm format
  strftime(buffer,80,"%Y-%m-%d_%H-%M",timeinfo);
  return std::string(buffer);
}

// reads in csv files
inline std::vector<std::vector<std::string>> CSVReaderString(std::string filename) {
  std::ifstream fin(filename);
  std::vector<std::vector<std::string>> data;

  std::string temp;

  while (fin >> temp) {

    std::vector<std::string> row;
    std::string line, word;

    std::getline(fin, line);
    std::istringstream s(line);

    while (std::getline(s, word, ',')) {
      row.push_back(word);
    }

    // add row to data.
    data.push_back(row);
  }

  return data;
}

// reads in csv files
template <typename T>
std::vector<std::vector<T>> CSVReader(std::string filename) {
  std::ifstream fin(filename);

  if (!std::is_same<T, double>::value && !std::is_same<T, float>::value && !std::is_same<T, int>::value) {
    T test;
    std::cout << "Specified type: " << typeid(test).name() << std::endl;
    throw std::invalid_argument("io_utils::CSVReader: Only supports type int, float, and double! Exiting.");
  }
  std::vector<std::vector<T>> data;


  while (fin) {

    std::vector<T> row;
    std::string line, word;

    std::getline(fin, line, '\n');
    if (line.empty())
      continue;

    std::istringstream s(line);

    while (std::getline(s, word, ',')) {
      // std::cout << word.c_str() << ", ";
      if (std::is_same<T, double>::value)
        row.push_back(std::stod(word));
      else if (std::is_same<T, float>::value)
        row.push_back(std::stof(word));
      else if (std::is_same<T, int>::value)
        row.push_back(std::stoi(word));
    }

    // add row to data.
    data.push_back(row);
  }

  return data;
}

// prints out parsed  csv files
template <typename T>
void printCSV(const std::vector<std::vector<T>>& data) {
  for (auto& line : data) {
    for (auto& word: line) {
      std::cout << word << ", ";
    }
    std::cout << "endl" << std::endl;
  }
  return;
}

} // namespace io
