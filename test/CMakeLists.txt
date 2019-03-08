include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../cpp ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
add_executable(test test.cc)
target_link_libraries(test ulb)