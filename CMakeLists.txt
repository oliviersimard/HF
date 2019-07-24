cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(HubbardTests)

# important for the use of tuples
SET(CMAKE_CXX_FLAGS "-std=c++0x")

set(YOUR_DIRECTORY "$ENV{HOME}/include")
set(YOUR_LIB_DIRECTORY "$ENV{HOME}/lib")

#Can manually add the sources using the set command as follows:
#set(SOURCES src/Hubbard.cpp)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB YOUR_LIBRARIES "${YOUR_LIB_DIRECTORY}/*json*.a")
file(GLOB SOURCES "src/*.cpp" "${YOUR_DIRECTORY}/*.h" ${YOUR_LIBRARIES})
#file(GLOB YOUR_LIBRARIES "${YOUR_LIB_DIRECTORY}/*json*.a")
set(LIBRARIES_CUSTOM fftw3 armadillo)

add_executable(runTestsHubbard.a ${SOURCES})

# Locate GTest
find_package(GTest REQUIRED)

if (GTest_FOUND)
    #Bring the headers into the project
    include_directories(${GTEST_INCLUDE_DIRS} src ${YOUR_DIRECTORY})
    # Link runTests with what we want to test and the GTest and pthread library
    message(${GTEST_LIBRARIES})
    target_link_libraries(runTestsHubbard.a ${GTEST_LIBRARIES} ${YOUR_LIBRARIES} pthread ${LIBRARIES_CUSTOM})
endif (GTest_FOUND) 
