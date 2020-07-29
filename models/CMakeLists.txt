include(FetchContent)
FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse
  GIT_TAG v2.1
  )
FetchContent_MakeAvailable(argparse)

add_executable(classical classical.cpp)
target_link_libraries(classical PRIVATE stochastic argparse::argparse)
install(TARGETS classical RUNTIME DESTINATION bin COMPONENT lmscpp)

add_executable(skewness skewness.cpp)
target_link_libraries(skewness PRIVATE stochastic argparse::argparse)
install(TARGETS skewness RUNTIME DESTINATION bin COMPONENT lmscpp)
