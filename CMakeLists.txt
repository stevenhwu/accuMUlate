cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(accuMUlate)

# guard against in-source builds

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()



set(CMAKE_CXX_FLAGS  "-std=c++11 -g -Wvla")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmakeModule/")

find_package( Boost COMPONENTS program_options REQUIRED )
find_package( Eigen3 3 REQUIRED )
find_package( Bamtools REQUIRED )
find_package( gsl REQUIRED )

FIND_PACKAGE(ZLIB)
IF(ZLIB_FOUND)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
ENDIF(ZLIB_FOUND)


set(LIBS ${LIBS} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${Bamtools_LIBRARIES})
set(INCLUDE_DIR ${INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src ${Boost_INCLUDE_DIR} ${Bamtools_INCLUDE_DIR}  ${EIGEN3_INCLUDE_DIR})


include_directories("${CMAKE_SOURCE_DIR}/src/third-party/bamtools/src")
include_directories(${INCLUDE_DIR})

ADD_SUBDIRECTORY(src)
