set( library_name HDK_PolyStokes )
add_library( ${library_name} SHARED
    HDK_PolyStokes.C
    HDK_PolyStokes.h
    HDK_PolyStokesSolver.cpp
    HDK_PolyStokesSolver.h
    HDK_PolyStokesSolver_Classifier.cpp
    HDK_PolyStokesSolver_ConstructMatrixBlocks.cpp
    HDK_PolyStokesSolver_Preconditioners.cpp
    HDK_PolyStokesSolver_AssembleBlocks.cpp
    HDK_PolyStokesSolver_AssembleSystem.cpp
)
target_include_directories( ${library_name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include )
# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name} HDK_PolyStokes_Lib Houdini Eigen3::Eigen)
# Configure several common target properties, such as its output directory.
houdini_configure_target( ${library_name} )
