set(prjname eagle_test)

set(all_headers )

file (GLOB dir_2_sources *.cpp *.c)
source_group("" FILES ${dir_2_sources})
set(all_sources ${dir_2_sources} )

add_executable(${prjname} ${all_headers} ${all_sources})

set_target_properties(${prjname} PROPERTIES DEBUG_OUTPUT_NAME ${prjname}${DEBUG_TARGET_SUFFIX})
set_target_properties(${prjname} PROPERTIES RELEASE_OUTPUT_NAME ${prjname})

include_directories(${PROJECTS_SRC_ROOT}/condor)

target_link_libraries(${prjname} debug math${DEBUG_TARGET_SUFFIX} optimized math)
target_link_libraries(${prjname} debug ang${DEBUG_TARGET_SUFFIX} optimized ang)
target_link_libraries(${prjname} debug env${DEBUG_TARGET_SUFFIX} optimized env)
target_link_libraries(${prjname} debug acft${DEBUG_TARGET_SUFFIX} optimized acft)
target_link_libraries(${prjname} debug eagle${DEBUG_TARGET_SUFFIX} optimized eagle)
target_link_libraries(${prjname} debug jail${DEBUG_TARGET_SUFFIX} optimized jail)
target_link_libraries(${prjname} ${Boost_LIBRARIES})

if (WIN32)
 STRING(TOUPPER eagle LIB)
 add_definitions(-D${LIB}_DLL)
endif(WIN32)
target_link_libraries(${prjname} ${Boost_LIBRARIES})

#=================================================================================================================
# Link Snippet for OpenSceneGraph
#=================================================================================================================
foreach(A_LIB ${OPENSCENEGRAPH_LIBS_RELEASE})
    target_link_libraries(${prjname} optimized ${A_LIB})
endforeach()
foreach(A_LIB ${OPENSCENEGRAPH_LIBS_DEBUG})
    target_link_libraries(${prjname} debug ${A_LIB})
endforeach()

#=================================================================================================================
# Link Snippet for OSGearth
#=================================================================================================================
foreach(A_LIB ${OSGEARTH_LITE_LIBS_RELEASE})
    target_link_libraries(${prjname} optimized ${A_LIB})
endforeach()
foreach(A_LIB ${OSGEARTH_LITE_LIBS_DEBUG})
    target_link_libraries(${prjname} debug ${A_LIB})
endforeach()

#=================================================================================================================
# Link Snippet for OpenCV
#=================================================================================================================
foreach(A_LIB ${OPENCV_LIBS_RELEASE}) 
	target_link_libraries(${prjname} optimized ${A_LIB})
endforeach()
foreach(A_LIB ${OPENCV_LIBS_DEBUG}) 
	target_link_libraries(${prjname} debug ${A_LIB})
endforeach()

