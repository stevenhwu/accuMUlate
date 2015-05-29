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


set(CMAKE_CXX_FLAGS_RelWithDebInfo "-std=c++11 -O3 -Wall -Wextra -pedantic -Wno-deprecated -Wno-deprecated-declarations -Wno-deprecated")
#set(CMAKE_CXX_FLAGS_RelWithDebInfo "-std=c++11 -O3 -w")
#set(CMAKE_CXX_FLAGS_RelWithDebInfo "${CMAKE_CXX_FLAGS_RelWithDebInfo} -Wno-deprecated-declarations -Wno-deprecated" )



set(CMAKE_CXX_FLAGS_DEBUG  "-std=c++11 -g  -O0 -Wvla -Wall  -Wextra -pedantic ")
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


#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wextra-semi ")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wconversion -Wsign-conversion ")

##Everything
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic") #forget about c++98

##turn off some for now, should add back to final code
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-deprecated-declarations -Wno-deprecated" )
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-unused-private-field -Wno-unused-variable -Wno-unused-private-field")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-unused-parameter -Wno-unused-function")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-unreachable-code")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-conversion -Wno-sign-conversion -Wno-missing-noreturn") #Way too many from Eigen/Boost
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-float-equal -Wno-sign-compare")


#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  \"
#"-Wno-deprecated-declarations -Wno-deprecated "\
#-Wno-float-equal -Wno-sign-compare")
#-Wno-unused-private-field -Wno-unused-variable -Wno-unused-private-field \
#-Wno-unused-parameter -Wno-unused-function \
#-Wno-unreachable-code \
#-Wno-conversion -Wno-sign-conversion -Wno-missing-noreturn \

##
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-shadow")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-padded")

##
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-global-constructors -Wno-exit-time-destructors")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-weak-vtables -Wno-missing-prototypes")



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

message("CMAKE_SYSTEM_PREFIX_PATH ${CMAKE_SYSTEM_PREFIX_PATH}")
message("CMAKE_SYSTEM_INCLUDE_PATH ${CMAKE_SYSTEM_INCLUDE_PATH}")
message("CMAKE_SYSTEM_LIBRARY_PATH ${CMAKE_SYSTEM_LIBRARY_PATH}")
message("CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH}")
message("CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH}")


message("CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}")
SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ~/lib ~/include ${CMAKE_INSTALL_PREFIX} ~)
message(${CMAKE_INSTALL_PREFIX})
message("CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}")


########## FIND_PACKAGE ##########

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmakeModule/")
find_package( Boost COMPONENTS program_options thread system REQUIRED )
find_package( Eigen3 3 REQUIRED )
find_package( Bamtools REQUIRED )
#find_package( gsl REQUIRED )

FIND_PACKAGE(ZLIB)
IF(ZLIB_FOUND)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
ENDIF(ZLIB_FOUND)


set(LIBS ${LIBS} ${Boost_LIBRARIES} ${Bamtools_LIBRARIES}) #  ${GSL_LIBRARIES}
set(INCLUDE_DIR ${INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${Bamtools_INCLUDE_DIR}  ${EIGEN3_INCLUDE_DIR})
message(${Bamtools_INCLUDE_DIR})
message(${Bamtools_LIBRARIES})

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(SYSTEM ${INCLUDE_DIR})

#include_directories(SYSTEM "/usr/include/")
#include_directories(SYSTEM "/usr/include/google")
##include_directories("/usr/include/sparsehash")
#include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/src/third-party/bamtools/src")
#include_directories("${CMAKE_SOURCE_DIR}/src/third-party/bamtools/src")
#include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/src/third-party/bamtools/src/api")
#include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/src/third-party/bamtools/src/utils")

message("==================== BUILD TYPE ${CMAKE_BUILD_TYPE} ====================")
ADD_SUBDIRECTORY(src)

########### Find/add GTest #################

find_package (GTest) #not setup properly yet
IF(GTEST_FOUND)

#1
#ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/include/googletest-trunk/")
#SET(GTEST_ROOT "${CMAKE_SOURCE_DIR}/include/googletest-trunk")

  find_package(Threads)
  message("${Threads_FOUND}")
  IF(THREADS_FOUND)
    SET(GTEST_THREAD_LIB "-lpthread")
    message("${ZLIB_INCLUDE_DIR}  ${ZLIB_LIBRARIES}  ${ZLIB_FOUND}")
    message("${GTEST_THREAD_LIB}A${CMAKE_THREAD_LIBS_INIT}B${Threads}====")
  ENDIF(THREADS_FOUND)

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

# Defines the following variables:
#
#   GTEST_FOUND - Found the Google Testing framework
#   GTEST_INCLUDE_DIRS - Include directories
#
# Also defines the library variables below as normal
# variables.  These contain debug/optimized keywords when
# a debugging library is found.
#
#   GTEST_BOTH_LIBRARIES - Both libgtest & libgtest-main
#   GTEST_LIBRARIES - libgtest
#   GTEST_MAIN_LIBRARIES - libgtest-main
#
ENDIF(GTEST_FOUND)
