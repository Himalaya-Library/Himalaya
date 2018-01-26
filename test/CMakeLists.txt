find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(-DBOOST_TEST_DYN_LINK)

function(add_himalaya_test name)
  add_executable(
    ${name}
    ${name}.cpp
    )
  target_link_libraries(
    ${name}
    Himalaya_static
    DSZ_static
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    )
  add_test(NAME ${name} COMMAND ${name})
endfunction()

add_himalaya_test(test_EFT_expressions)
add_himalaya_test(test_lambda_limits)
add_himalaya_test(test_lambda_normalization)
