cmake_minimum_required(VERSION 3.4)
project (test)

enable_testing()
find_package(SymEngine 0.1.0 CONFIG PATHS /home/thiago/.local/lib/symengine/lib/cmake/symengine)
# find_package(SymEngine 0.1.0 CONFIG PATHS /home/thiago/github/symengine/build)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_compile_options(-Wall -Wextra -Winline -O0 -g3)


include_directories(${CMAKE_SOURCE_DIR}/lib/bitsery/include)
include_directories(${SYMENGINE_INCLUDE_DIRS})
add_executable(testsym sym.cpp)
target_link_libraries(testsym ${SYMENGINE_LIBRARIES})

add_executable(testbitsery bit.cpp)
# add_executable(testbitsery2 bit2.cpp)
add_executable(testbitsery3 bit3.cpp)
add_executable(lms main.cpp stochastic.cpp utils.cpp)
target_link_libraries(lms ${SYMENGINE_LIBRARIES})

add_executable(testall testall.cpp stochastic.cpp utils.cpp)
target_link_libraries(testall ${SYMENGINE_LIBRARIES})
add_test(testall testall)
