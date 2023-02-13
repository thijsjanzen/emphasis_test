#==============================================================================
# emphasis CMake file
#
# Hanno Hildenbrandt 2023
#============================================================================== 
cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED)
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")

if (WIN32)
    add_definitions(/DNOMINMAX /DWIN32_LEAN_AND_MEAN /D_CRT_SECURE_NO_WARNINGS ${GLSL_VERSION})
else()
    # ...
endif()

project("emphasis_test" CXX)

find_package(TBB CONFIG REQUIRED)

file(GLOB src
    ${PROJECT_SOURCE_DIR}/*.h
    ${PROJECT_SOURCE_DIR}/*.hpp
    ${PROJECT_SOURCE_DIR}/*.cpp
)

add_executable(emphasis ${src})
target_link_libraries(emphasis PRIVATE TBB::tbb)
set_target_properties(emphasis PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/$<0:>)

