add_definitions(-DDOCTEST_CONFIG_IMPLEMENT_WITH_MAIN)
add_definitions(-DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

function(add_himalaya_test name)
  add_executable(
    ${name}
    ${name}.cpp
    )
  target_link_libraries(
    ${name}
    Himalaya::Himalaya
    Himalaya::DSZ
    )
  target_include_directories(
    ${name} PRIVATE $<TARGET_PROPERTY:Himalaya::Himalaya,INCLUDE_DIRECTORIES>
    )
  add_test(NAME ${name} COMMAND ${name})
endfunction()

add_himalaya_test(test_benchmark)
add_himalaya_test(test_example)
add_himalaya_test(test_EFT_expressions)
add_himalaya_test(test_FO_expressions)
add_himalaya_test(test_FO_vs_EFT_expressions)
add_himalaya_test(test_Himalaya_interface)
add_himalaya_test(test_Numerics)
add_himalaya_test(test_li2)
add_himalaya_test(test_delta_lambda)
add_himalaya_test(test_points)
add_himalaya_test(test_PV)

if(Mathematica_WolframLibrary_FOUND)
  Mathematica_WolframLibrary_ADD_TEST (
    NAME test_LibraryLink
    TARGET Himalaya_LibraryLink
    SCRIPT test_LibraryLink.m)
endif()
