#==============================================================================
# emphasis CMake file
#
# Hanno Hildenbrandt 2020
#============================================================================== 
cmake_minimum_required(VERSION 3.20)

project(emphasis_test)

file(GLOB src
    ${PROJECT_SOURCE_DIR}/*.h
    ${PROJECT_SOURCE_DIR}/*.hpp
    ${PROJECT_SOURCE_DIR}/*.cpp
)

add_executable(emphasis ${src})
set_target_properties(emphasis PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/$<0:>)

