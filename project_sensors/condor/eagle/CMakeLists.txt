set(prjname eagle)

file (GLOB dir_0_headers *.h *.hpp)
SOURCE_GROUP("" FILES ${dir_0_headers})
set(all_headers ${dir_0_headers} )

file (GLOB dir_1_sources *.cpp *.c)
SOURCE_GROUP("" FILES ${dir_1_sources})
set(all_sources ${dir_1_sources} )

STRING(TOUPPER ${prjname} PRJNAME)
add_definitions(-D${PRJNAME}_DLL)
add_definitions(-D${PRJNAME}_DLL_EXPORTS)

add_library(${prjname} SHARED ${all_headers} ${all_sources})

set_target_properties(${prjname} PROPERTIES DEBUG_OUTPUT_NAME ${prjname}${DEBUG_TARGET_SUFFIX})
set_target_properties(${prjname} PROPERTIES RELEASE_OUTPUT_NAME ${prjname})

include_directories(${PROJECTS_SRC_ROOT}/condor)

target_link_libraries(${prjname} debug ang${DEBUG_TARGET_SUFFIX} optimized ang)
target_link_libraries(${prjname} debug env${DEBUG_TARGET_SUFFIX} optimized env)
target_link_libraries(${prjname} debug acft${DEBUG_TARGET_SUFFIX} optimized acft)
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


