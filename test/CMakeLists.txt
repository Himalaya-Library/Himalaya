find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
add_definitions (-DBOOST_TEST_DYN_LINK)

add_executable(test_lambda_normalization test_lambda_normalization.cpp)
target_link_libraries(
  test_lambda_normalization
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(NAME test_lambda_normalization COMMAND test_lambda_normalization)
