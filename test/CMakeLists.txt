find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(-DBOOST_TEST_DYN_LINK)

# test degenerate limits for zeta_lambda
add_executable(
  test_lambda_limits
  test_lambda_limits.cpp
)
target_link_libraries(
  test_lambda_limits
  Himalaya_static
  DSZ_static
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(NAME test_lambda_limits COMMAND test_lambda_limits)
