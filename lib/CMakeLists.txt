include(FetchContent)
FetchContent_Declare(
  symengine
  GIT_REPOSITORY https://github.com/symengine/symengine
  # GIT_TAG v0.6.0
  )
FetchContent_MakeAvailable(symengine)
find_package(SymEngine)

add_library(stochastic stochastic.cpp utils.cpp )
target_include_directories(stochastic PUBLIC ../include)
target_include_directories(stochastic PUBLIC ${SYMENGINE_INCLUDE_DIRS})
target_link_libraries(stochastic PUBLIC ${SYMENGINE_LIBRARIES})
