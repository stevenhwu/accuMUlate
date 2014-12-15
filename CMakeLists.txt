cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(accuMUlate)

# guard against in-source builds

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)


set(CMAKE_CXX_FLAGS  "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE  "-std=c++11 -O3")
set(CMAKE_CXX_FLAGS_RelWithDebInfo "-std=c++11 -O3 -Wall")
set(CMAKE_CXX_FLAGS_DEBUG  "-std=c++11 -g -O0 -Wvla -Wall  -Wextra -pedantic ")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} "-std=c++11 -g -Wvla -Wall  -Wextra -pedantic -Wno-unused-private-field -Wno-unused-variable")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wextra -Waggregate-return -Wcast-align") #-pedantic-errors
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wcast-qual  -Wchar-subscripts  -Wcomment " )
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wdisabled-optimization ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wfloat-equal  -Wformat  -Wformat=2 ")#-Werror
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wformat-nonliteral -Wformat-security  ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wformat-y2k ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wimplicit  -Wimport  -Winit-self  -Winline ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Winvalid-pch   ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wlong-long -Wmissing-braces ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wmissing-field-initializers -Wmissing-format-attribute   ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wmissing-include-dirs ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wpacked -Wparentheses  -Wpointer-arith ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wredundant-decls -Wreturn-type ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wsequence-point   -Wsign-compare  -Wstack-protector ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wswitch-default ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wswitch-enum -Wtrigraphs  -Wuninitialized ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wunknown-pragmas  -Wunused -Wunreachable-code")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wunused-function  -Wunused-label  -Wunused-parameter ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wunused-value  -Wunused-variable  -Wvariadic-macros ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wvolatile-register-var  -Wwrite-strings ")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-unused-private-field -Wno-unused-variable -Wno-unused-private-field")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-unreachable-code") #for debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-conversion -Wno-padded -Wno-missing-noreturn") #Way too many from Eigen/Boost
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wshadow")


#set(CMAKE_CXX_FLAGS_DEBUG  "-std=c++11 -g -Wvla -Wall  -Wextra ")



#set(CMAKE_CXX_FLAGS  "-std=c++11 -g -Wvla -Weverything -Wno-unused-private-field -Wno-unused-variable -Wno-documentation -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-deprecated -Wno-sign-conversion")

########## BUILD TYPE ##########
SET( CMAKE_CXX_FLAGS_MAINTAINER "-Wall -Wabi" CACHE STRING
    "Flags used by the C++ compiler during maintainer builds."
    FORCE )
SET( CMAKE_C_FLAGS_MAINTAINER "-Wall -pedantic" CACHE STRING
    "Flags used by the C compiler during maintainer builds."
    FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_MAINTAINER
    "-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
    "Flags used for linking binaries during maintainer builds."
    FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
    "-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
    "Flags used by the shared libraries linker during maintainer builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_MAINTAINER
    CMAKE_C_FLAGS_MAINTAINER
    CMAKE_EXE_LINKER_FLAGS_MAINTAINER
    CMAKE_SHARED_LINKER_FLAGS_MAINTAINER )
# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Maintainer."
    FORCE )




########## FIND_PACKAGE ##########

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

include_directories("/usr/include/")
include_directories("/usr/include/google")
include_directories("${CMAKE_SOURCE_DIR}/src/third-party/bamtools/src")
include_directories(${INCLUDE_DIR})


message("------------------------BUILD TYPE ${CMAKE_BUILD_TYPE}")
ADD_SUBDIRECTORY(src)

########### Find/add GTest #################
## find_package (GTEST) not setup properly yet
#1
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/include/googletest-trunk/")
#SET(GTEST_ROOT "${CMAKE_SOURCE_DIR}/include/googletest-trunk")
ADD_SUBDIRECTORY(test)


#2 Some link error??
#SET(GTEST_ROOT "${CMAKE_SOURCE_DIR}/include/googletest-trunk")
#find_package( GTest  )
#IF(GTEST_FOUND)
#    message("Found GTEST at ${GTEST_INCLUDE_DIRS} ${gtest_SOURCE_DIR} ")
#    include_directories(${GTEST_INCLUDE_DIRS})
#    ADD_SUBDIRECTORY(test)
#ENDIF(GTEST_FOUND)

#3 SVN checkout version
#ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/gtest")







