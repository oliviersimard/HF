cmake_minimum_required(VERSION 3.02 FATAL_ERROR)

project(HubbardTests)

# important for the use of tuples
SET(CMAKE_CXX_FLAGS "-std=c++0x")

set(LIBRARIES_CUSTOM fftw3 armadillo json_spirit)

set(SOURCES_COMMON
        src/fft.cpp src/Hubbard.cpp
        src/HubbardSelfCon.cpp
        src/json_utils.cpp src/param.cpp)

set(SOURCES_HF hartree2DSure.cpp ${SOURCES_COMMON})
set(SOURCES_MAIN src/Main.cpp ${SOURCES_COMMON})

#------Executables
add_executable(HF ${SOURCES_HF})
target_link_libraries(HF ${LIBRARIES_CUSTOM})
add_executable(MAIN ${SOURCES_MAIN})
target_link_libraries(MAIN ${LIBRARIES_CUSTOM})


# Taken from github google test repo

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif()

# Now simply link against gtest or gtest_main as needed. Eg
enable_testing()
add_executable(hubbard_test test/Hubbard_unit_tests.cpp ${SOURCES_COMMON})
target_link_libraries(hubbard_test ${LIBRARIES_CUSTOM} gtest_main gtest)
add_test(NAME hubbard_test COMMAND ./hubbard_test)
