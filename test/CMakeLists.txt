include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../cpp)
add_executable(test test.cc)
target_link_libraries(test ulb)