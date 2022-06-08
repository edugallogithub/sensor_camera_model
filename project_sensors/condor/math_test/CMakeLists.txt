
set(prjname math_test)

file (GLOB dir_27_headers *.h *.hpp)
source_group("" FILES ${dir_27_headers})
set(all_headers ${dir_27_headers} )

file (GLOB dir_28_sources *.cpp *.c)
source_group("" FILES ${dir_28_sources})
set(all_sources ${dir_28_sources} )

add_executable(${prjname} ${all_headers} ${all_sources})

set_target_properties(${prjname} PROPERTIES DEBUG_OUTPUT_NAME ${prjname}${DEBUG_TARGET_SUFFIX})
set_target_properties(${prjname} PROPERTIES RELEASE_OUTPUT_NAME ${prjname})

include_directories(${PROJECTS_SRC_ROOT}/condor)

target_link_libraries(${prjname} debug math${DEBUG_TARGET_SUFFIX} optimized math)
target_link_libraries(${prjname} debug jail${DEBUG_TARGET_SUFFIX} optimized jail)
target_link_libraries(${prjname} ${Boost_LIBRARIES})

foreach(A_LIB ${EIGEN_LIBS_RELEASE}) 
	target_link_libraries(${prjname} optimized ${A_LIB}) # adds all Eigen libraries to target
endforeach()
foreach(A_LIB ${EIGEN_LIBS_DEBUG}) 
	target_link_libraries(${prjname} debug ${A_LIB}) # adds all Eigen libraries to target
endforeach()
