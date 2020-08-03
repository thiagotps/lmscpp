cmake_minimum_required(VERSION 3.16)

project (lmscpp VERSION 0.1
  DESCRIPTION "An implementation of Douglas Scott's exact expectation analysis"
  LANGUAGES CXX)

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-march=native -mtune=native")
set(CMAKE_C_FLAGS "-march=native -mtune=native")
set (CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set (CMAKE_C_FLAGS_RELEASE "-Ofast")

enable_testing()

add_subdirectory(lib)
add_subdirectory(tests)
add_subdirectory(models)
add_subdirectory(samples)
add_subdirectory(scripts)
