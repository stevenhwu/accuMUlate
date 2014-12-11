


SET(SOURCE_FILES CACHE INTERNAL "Global")
message("A: ${SOURCE_FILES}")
SET(SOURCE_FILES ${SOURCE_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/model.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/parsers.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/VariantVisitor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/variant_visitor_two.cc

    ${CMAKE_CURRENT_SOURCE_DIR}/lookup.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/mutation_prob.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/site_prob.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/sequence_prob.cc


    CACHE INTERNAL "Global"
    )

SET(SOURCE_FILES_accuMUlate ${SOURCE_FILES_accuMUlate}
    ${CMAKE_CURRENT_SOURCE_DIR}/model.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/parsers.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/VariantVisitor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/variant_visitor_two.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lookup.cc
    )

#list (APPEND SOURCE_LIST ${SOURCE_FILES})
#file(GLOB_RECURSE SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

message("BB: ${SOURCE_FILES} ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ")

ADD_SUBDIRECTORY(distributions)
ADD_SUBDIRECTORY(evolution_models)
ADD_SUBDIRECTORY(algorithm)


message("C: ${SOURCE_FILES} ")

add_executable(accuMUlate main.cc ${SOURCE_FILES_accuMUlate})
target_link_libraries(accuMUlate ${LIBS} )

add_executable(EMTest EMTest.cc ${SOURCE_FILES})
target_link_libraries(EMTest ${LIBS} )
target_link_libraries(EMTest profiler)
#target_link_libraries(accuMUlate distributions)
#target_link_libraries(accuMUlate bamtools bamtools-utils ${LIBS} distributions)


#add_executable(pp utils/post_processor.cc parsers.cc model.cc)
#target_link_libraries(pp ${Boost_LIBRARIES} ${Bamtools_LIBRARIES})
#link_directories("${Bamtools_LIBRARY_DIR}/bamtools/")
#target_link_libraries(pp ${Boost_LIBRARIES} bamtools bamtools-utils )


#add_executable(freq_seq utils/freq_seq.cc parsers.cc)
#target_link_libraries(pp bamtools bamtools-utils)
