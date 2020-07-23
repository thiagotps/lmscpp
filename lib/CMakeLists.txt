# Prevent SymEngine from building tests
set (BUILD_TESTS CACHE BOOL "NO")

include(FetchContent)
FetchContent_Declare(
  symengine
  GIT_REPOSITORY https://github.com/symengine/symengine
  # GIT_TAG v0.6.0
  )
FetchContent_MakeAvailable(symengine)

FetchContent_Declare(
  bitsery
  GIT_REPOSITORY https://github.com/fraillt/bitsery
  GIT_TAG v5.1.0
  )

FetchContent_MakeAvailable(bitsery)

find_package(SymEngine)
find_package(Threads)

add_library(stochastic stochastic.cpp utils.cpp)
target_compile_features(stochastic PUBLIC cxx_std_17)
target_include_directories(stochastic PUBLIC ../include ${SYMENGINE_INCLUDE_DIRS})
target_link_libraries(stochastic PUBLIC ${SYMENGINE_LIBRARIES} Bitsery::bitsery ${CMAKE_THREAD_LIBS_INIT})
