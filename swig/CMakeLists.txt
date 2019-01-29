#include(SwigPython)
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(NumPy REQUIRED)
include_directories(${NUMPY_INCLUDE_DIRS})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cpp)

set_source_files_properties(ulb.i PROPERTIES CPLUSPLUS ON)
swig_add_library(pyulb TYPE SHARED LANGUAGE python SOURCES ulb.i)
swig_link_libraries(pyulb ${PYTHON_LIBRARIES} ulb)
#add_swig_python_module(pyulb ulb.i
#        CPLUSPLUS
#        INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../cpp
#        LINK_LIBRARIES ulb
#        SWIG_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
#        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pyulb)
