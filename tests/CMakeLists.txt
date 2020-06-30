add_executable(testall testall.cpp)
target_link_libraries(testall PRIVATE stochastic)
add_test(NAME testall COMMAND testall)
