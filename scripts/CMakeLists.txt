add_test(NAME classical COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/classicaltest.sh $<TARGET_FILE:classical> ${CMAKE_SOURCE_DIR}/scripts/base_results)
add_test(NAME skewness COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/skewness_test.sh $<TARGET_FILE:skewness> ${CMAKE_SOURCE_DIR}/scripts/skewness_results)
