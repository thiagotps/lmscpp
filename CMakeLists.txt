cmake_minimum_required(VERSION 3.17)

project (lmscpp VERSION 0.1
  DESCRIPTION "An implementation of Douglas Scott's exact expectation analysis"
  LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()

add_subdirectory(lib)
add_subdirectory(tests)
