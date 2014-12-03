
################################
# GTest
################################

enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

#/home/steven/Postdoc2/Project_accuMUlate/accuMUlate/include/gtest-1.7.0

################################
# Unit Tests
################################
# Add test cpp file
#SET(TEST_SOURCE_FILES ${TEST_SOURCE_FILES} testA.cc testB.cc)
SET(TEST_SOURCE_FILES ${TEST_SOURCE_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/MutationProbTest.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/SequenceProbTest.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/SiteProbTest.cc
    )

ADD_SUBDIRECTORY("evolution_models")

add_executable( UnitTests ${TEST_SOURCE_FILES} ${SOURCE_FILES})


# Link test executable against gtest & gtest_main
target_link_libraries(UnitTests gtest gtest_main ${LIBS} )
add_test( UnitTests UnitTests )


