add_executable(testall testall.cpp)
target_link_libraries(testall PRIVATE stochastic)
add_test(NAME testall COMMAND testall)

add_executable(testnode testnode.cpp)
target_link_libraries(testnode PRIVATE stochastic)
add_test(NAME testnode COMMAND testnode)

add_executable(testmatrix testmatrix.cpp)
target_link_libraries(testmatrix PRIVATE stochastic)
add_test(NAME testmatrix COMMAND testmatrix)
