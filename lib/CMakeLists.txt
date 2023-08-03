set( library_name HDK_PolyStokes_Lib )
set( library_include_dir ${CMAKE_CURRENT_SOURCE_DIR}/include )
set( library_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_library( ${library_name} STATIC
	${library_source_dir}/ApplyPressureStressMatrix.cpp
	${library_source_dir}/Preconditioner.cpp
	#${library_source_dir}/GridClassifier.cpp
)
#set_target_properties( ${library_name} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories( ${library_name} PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include )

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name} Houdini Eigen3::Eigen)