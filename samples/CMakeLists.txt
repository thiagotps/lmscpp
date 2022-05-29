include(FetchContent)
FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse
  GIT_TAG v2.5
  )
FetchContent_MakeAvailable(argparse)

FetchContent_Declare(
  pcg
  GIT_REPOSITORY https://github.com/imneme/pcg-cpp
  GIT_TAG v0.98.1
  )
FetchContent_MakeAvailable(pcg)

add_executable(skempold skempold.cpp)
target_link_libraries(skempold PRIVATE stochastic argparse::argparse)
install(TARGETS skempold RUNTIME DESTINATION bin COMPONENT lmscpp)

add_executable(skemp skemp.cpp)
target_include_directories(skemp PRIVATE "${pcg_SOURCE_DIR}/include")
target_link_libraries(skemp PRIVATE stochastic argparse::argparse)
install(TARGETS skemp RUNTIME DESTINATION bin COMPONENT lmscpp)

add_executable(fourthemp fourthemp.cpp)
target_include_directories(fourthemp PRIVATE "${pcg_SOURCE_DIR}/include")
target_link_libraries(fourthemp PRIVATE stochastic argparse::argparse)
install(TARGETS fourthemp RUNTIME DESTINATION bin COMPONENT lmscpp)
