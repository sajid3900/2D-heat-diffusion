cmake_minimum_required(VERSION 3.5)

project(SiSc_lab2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#add_executable(SiSc_lab2d SiSc_lab2d.cpp )

INCLUDE_DIRECTORIES(dcw6i34nel/include/  eigen-3.3.7/)
LINK_DIRECTORIES(dcw6i34nel/lib)
add_executable(SiSc_lab2d main.cpp)
if(WIN32)
TARGET_LINK_LIBRARIES(SiSc_lab2d libdcoc_MD.lib)
else()
TARGET_LINK_LIBRARIES(SiSc_lab2d libdcoc.a)
endif()
#find_library(UMLPACK)
#target_link_libraries(SiSc_lab2d myfancyLib)

#set(SOURCEFILES SiSc_lab2d.cpp test2.cpp)
#add_executable(SiSc_lab2d_123 ${SOURCEFILES})

include_directories(.)

############ GTEST
# Download and unpack googletest at configure time
configure_file(${CMAKE_SOURCE_DIR}/additionals/gtest.cmake googletest-download/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
    ${CMAKE_BINARY_DIR}/googletest-build
    EXCLUDE_FROM_ALL)

########## END GTEST


include(${CMAKE_SOURCE_DIR}/additionals/cmake_scripts/stce_utils.cmake)


#add tests target for recursively adding subfolders
stce_add_recursively(${CMAKE_SOURCE_DIR}/tests/)


