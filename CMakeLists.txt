cmake_minimum_required(VERSION 3.16)

project (lmscpp VERSION 0.1
  DESCRIPTION "An implementation of Douglas Scott's exact expectation analysis"
  LANGUAGES CXX)

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()

add_subdirectory(lib)
add_subdirectory(tests)
add_subdirectory(models)
add_subdirectory(samples)
add_subdirectory(scripts)
add_subdirectory(utils)
add_subdirectory(deps/eigen-3.4-rc1)
add_subdirectory(deps/spectra-1.0.0)
