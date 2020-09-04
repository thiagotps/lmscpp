include(FetchContent)

FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse
  GIT_TAG v2.1
  )
FetchContent_MakeAvailable(argparse)

add_executable(gauss_emp gauss_emp.cpp)
target_link_libraries(gauss_emp PRIVATE argparse::argparse)
target_compile_features(gauss_emp PUBLIC cxx_std_17)
install(TARGETS gauss_emp RUNTIME DESTINATION bin COMPONENT lmscpp)
