include(FetchContent)
FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse
  GIT_TAG v2.1
  )
FetchContent_MakeAvailable(argparse)

add_executable(skemp skemp.cpp)
target_link_libraries(skemp PRIVATE stochastic argparse::argparse)
install(TARGETS skemp RUNTIME DESTINATION bin COMPONENT lmscpp)

add_executable(skemp2 skemp2.cpp)
target_link_libraries(skemp2 PRIVATE stochastic argparse::argparse)
install(TARGETS skemp2 RUNTIME DESTINATION bin COMPONENT lmscpp)
