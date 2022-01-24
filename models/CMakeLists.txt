include(FetchContent)
FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse
  GIT_TAG v2.1
  )
FetchContent_MakeAvailable(argparse)

add_executable(classical classical.cpp common.cpp)
target_link_libraries(classical PRIVATE stochastic argparse::argparse)
install(TARGETS classical RUNTIME DESTINATION bin COMPONENT lmscpp)

add_executable(skewness skewness.cpp common.cpp)
target_link_libraries(skewness PRIVATE stochastic argparse::argparse)
install(TARGETS skewness RUNTIME DESTINATION bin COMPONENT lmscpp)

add_executable(fourth-order fourth-order.cpp common.cpp)
target_link_libraries(fourth-order PRIVATE stochastic argparse::argparse)
target_include_directories(fourth-order PRIVATE ../deps/eigen-3.4.0 ../deps/spectra-1.0.0/include)
install(TARGETS fourth-order RUNTIME DESTINATION bin COMPONENT lmscpp)
