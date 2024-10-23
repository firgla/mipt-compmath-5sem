include(FetchContent)
cmake_minimum_required(VERSION 3.28)
set(CMAKE_CXX_STANDARD 17)



set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 REQUIRED)

include_directories(EIGEN3_INCLUDE_DIR)

message ("${EIGEN3_INCLUDE_DIR}")

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.11.0
)
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)


file(GLOB_RECURSE TEST_SOURCES test*.cpp)

foreach (TEST_FILE IN LISTS TEST_SOURCES)
    get_filename_component(TEST_FILE_NAME ${TEST_FILE} NAME_WE)
    set(EXECUTABLE_NAME run_${TEST_FILE_NAME})
    add_executable(${EXECUTABLE_NAME} ${TEST_FILE})
    target_compile_options(${EXECUTABLE_NAME} PRIVATE -Werror -Wall -Wextra -pedantic -Wshadow -Wconversion -Wsign-conversion)
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${PROJECT_NAME} GTest::gtest_main)
    add_test(NAME ${EXECUTABLE_NAME} COMMAND ${EXECUTABLE_NAME})
endforeach ()