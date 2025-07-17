cmake_minimum_required(VERSION 3.31)
project(mpm_2d)

set(CMAKE_CXX_STANDARD 20)

add_executable(mpm_2d src/main.cpp)

# --- Find Packages ---
find_package(Eigen3 REQUIRED)
find_package(HDF5 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

# --- Collect Source Files ---
file(GLOB_RECURSE SOURCES "src/*.cpp")

# --- Create Executable ---
add_executable(${PROJECT_NAME} ${SOURCES})

# --- Configure Target ---
# Add include directories
target_include_directories(${PROJECT_NAME} PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
        ${EIGEN3_INCLUDE_DIR}
        ${HDF5_INCLUDE_DIRS}
)

# Link all necessary libraries
target_link_libraries(${PROJECT_NAME} PUBLIC
        ${HDF5_LIBRARIES}
        OpenGL::GL
        GLUT::GLUT
)