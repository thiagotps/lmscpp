FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse
  GIT_TAG v2.1
  )
FetchContent_MakeAvailable(argparse)

add_executable(classical classical.cpp)
target_link_libraries(classical PRIVATE stochastic argparse::argparse)
# add_test(NAME testnode COMMAND testnode)
