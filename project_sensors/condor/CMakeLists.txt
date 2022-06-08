
set(prjname condor_member)
set(version 1)

project(${prjname} CXX)
add_subdirectory(jail)
add_subdirectory(math)
add_subdirectory(math_test)
add_subdirectory(ang)
add_subdirectory(ang_test)
add_subdirectory(env)
add_subdirectory(env_test)
add_subdirectory(acft)
add_subdirectory(acft_test)
add_subdirectory(eagle)
add_subdirectory(eagle_test)


