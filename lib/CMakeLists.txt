# Prevent SymEngine from building tests
set (BUILD_TESTS OFF)

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

add_library(stochastic stochastic.cpp utils.cpp )
target_include_directories(stochastic PUBLIC ../include)
target_include_directories(stochastic PUBLIC ${SYMENGINE_INCLUDE_DIRS})
target_link_libraries(stochastic PUBLIC ${SYMENGINE_LIBRARIES})
target_link_libraries(stochastic PUBLIC Bitsery::bitsery)
