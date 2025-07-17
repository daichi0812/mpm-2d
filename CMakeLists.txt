cmake_minimum_required(VERSION 3.31)
project(mpm_2d LANGUAGES CXX)

# Require
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type if none supplied
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

#
# Options
#
option(MPM2D_BUILD_EXAMPLES "Build example / test programs" OFF)
option(MPM2D_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

#
# Dependencies
#

# Eigen3 (header-only)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# HDF5
find_package(HDF5 REQUIRED COMPONENTS C HL)

# OpenGL & GLUT (Rendering)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

# RapidXML (header-only)
# Try common locations + project external/rapidxml
find_path(RAPIDXML_INCLUDE_DIR
        NAMES rapidxml.hpp
        HINTS
        ${CMAKE_CURRENT_SOURCE_DIR}/external/rapidxml
        /opt/homebrew/include
        /usr/local/include
        /usr/include
        PATH_SUFFIXES rapidxml
)
if(NOT RAPIDXML_INCLUDE_DIR)
    message(FATAL_ERROR "rapidxml.hpp not found. Please download RapidXML and place it in external/rapidxml or install system-wide.")
endif()

#
# Sources
#
file(GLOB_RECURSE MPM2D_SOURCES
        CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
)

file(GLOB_RECURSE MPM2D_HEADERS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
)

add_executable(${PROJECT_NAME}
        ${MPM2D_SOURCES}
        ${MPM2D_HEADERS}   # not strictly needed, but lets IDEs show headers
)

#
# Include directories
#
target_include_directories(${PROJECT_NAME}
        PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
        ${EIGEN3_INCLUDE_DIRS}    # Provided by FindEigen3
        ${HDF5_INCLUDE_DIRS}
        ${RAPIDXML_INCLUDE_DIR}
)

#
# Compile options
#
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
    if(MPM2D_WARNINGS_AS_ERRORS)
        target_compile_options(${PROJECT_NAME} PRIVATE /WX)
    endif()
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
    if(MPM2D_WARNINGS_AS_ERRORS)
        target_compile_options(${PROJECT_NAME} PRIVATE -Werror)
    endif()
endif()

# macOS: ソース側で GL_SILENCE_DEPRECATION を定義しているので
# ここでは追加定義しない（警告回避）。必要ならコメントアウトを外す。
# if(APPLE)
#     target_compile_definitions(${PROJECT_NAME} PRIVATE GL_SILENCE_DEPRECATION)
# endif()

# Link libraries
target_link_libraries(${PROJECT_NAME}
        PUBLIC
        Eigen3::Eigen         # header-only, but target OK
        ${HDF5_LIBRARIES}
        OpenGL::GL
        GLUT::GLUT
)


# Some HDF5 builds require extra vars
if(HDF5_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${HDF5_C_LIBRARY_hdf5})
endif()

# RPATH convenience on *nix so you can run from build tree if any shared libs
if(UNIX AND NOT APPLE)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            BUILD_RPATH "$ORIGIN"
            INSTALL_RPATH "$ORIGIN"
    )
endif()

# ------------------------------------------------------------------
# Install (optional)
# ------------------------------------------------------------------
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)