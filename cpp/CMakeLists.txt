cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(sgt LANGUAGES CXX)

find_package(OpenMP)

set (RUN_WATERBOX_SRC
    ./IO/readInputFiles.cpp
    ./Quadrature/hermiteRootsWeights.cpp
    ./Quadrature/sphericalBessel.cpp
    ./Atom/utility.cpp
    ./Atom/atom.cpp
    ./ErfPointCharge/erfPointCharge.cpp
    ./Cube/cube.cpp
    ./Cube/vecOfCubes.cpp
    ./Hierarchy/hierarchy.cpp
    ./TestCode/run_waterBox.cpp
)

add_executable(run_waterbox ${RUN_WATERBOX_SRC})

set (DOUBLE_ATOMS_SRC
    ./IO/readInputFiles.cpp
    ./Atom/utility.cpp
    ./Atom/atom.cpp
    ./TestCode/doubleAtomsInX_Y_or_Z.cpp
)

add_executable(double_waterbox ${DOUBLE_ATOMS_SRC})

set_target_properties(
    run_waterbox
    double_waterbox
    PROPERTIES CXX_STANDARD 17
)

target_link_libraries(
    run_waterbox
    OpenMP::OpenMP_CXX
)

# Make sure that Boost can be found
include_directories(./Quadrature)
