# compile using "mkdir build && cd build && cmake .. -G "Visual Studio 16 2019"
# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.9 )
set (CMAKE_CXX_STANDARD 17)
project( HDK_PolyStokes )

# CMAKE_PREFIX_PATH must contain the path to the toolkit/cmake subdirectory of 
# the Houdini installation. See the "Compiling with CMAKE" section of the HDK
# documentation for more details, which describes several options
# for specifying this path.
list( APPEND CMAKE_PREFIX_PATH "$ENV{HFS}/toolkit/cmake" )

# Locate Houdini's libraries and header files.
# Registers an imported library target named 'Houdini'.
find_package( Houdini REQUIRED )

set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/extern/eigen")
add_subdirectory(${EIGEN3_INCLUDE_DIR})
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(lib)
add_subdirectory(exec)