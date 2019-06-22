find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(NumPy REQUIRED)
include_directories(${NUMPY_INCLUDE_DIRS})

find_package(PythonInterp REQUIRED)

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cpp)

set(genmeshattrs_py_in "${CMAKE_CURRENT_SOURCE_DIR}/genmeshattrs.py.in")
set(genmeshattrs_py    "${CMAKE_CURRENT_SOURCE_DIR}/genmeshattrs.py")
set(meshattrs_in "${CMAKE_CURRENT_SOURCE_DIR}/meshattr.i")
set(meshattrs "${CMAKE_CURRENT_SOURCE_DIR}/meshattrs.yaml")
configure_file(${genmeshattrs_py_in} ${genmeshattrs_py})

set_source_files_properties(ulb.i PROPERTIES CPLUSPLUS ON)
swig_add_library(pyulb TYPE SHARED LANGUAGE python SOURCES ulb.i)
swig_link_libraries(pyulb ${PYTHON_LIBRARIES} ulb)
#add_dependencies(pyulb  ${CMAKE_CURRENT_SOURCE_DIR}/meshattrs.yaml)

add_custom_command(TARGET pyulb
        PRE_BUILD
        #DEPENDS
        #meshattrs.yaml
        COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/genmeshattrs.py
        )
