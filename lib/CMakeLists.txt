# Prevent SymEngine from building tests
set (BUILD_TESTS CACHE BOOL "NO")

include(FetchContent)
FetchContent_Declare(
  symengine
  GIT_REPOSITORY https://github.com/symengine/symengine
  # NOTE: There is something preventing version v0.6.0 of being properly configured by CMAKE
  # Menwhile we will just use the master version and wait until a new stable version is released
  GIT_TAG 2ed68bc890a2cbe8507fec1b0a2a833d75b00dbe
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
