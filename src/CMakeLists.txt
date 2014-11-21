
message("A: ${SOURCE_FILES}")

SET(SOURCE_FILES  model.cc parsers.cc SequenceProb.cc MutationProb.cc VariantVisitor.cc Lookup.cc)


ADD_SUBDIRECTORY(distributions)
ADD_SUBDIRECTORY(models)

message("B: ${SOURCE_FILES}")


message("C: ${SOURCE_FILES}")

add_executable(accuMUlate main.cc ${SOURCE_FILES})
target_link_libraries(accuMUlate ${LIBS} )

add_executable(EMTest EMTest.cc ${SOURCE_FILES})
target_link_libraries(EMTest ${LIBS} )

#target_link_libraries(accuMUlate distributions)
#target_link_libraries(accuMUlate bamtools bamtools-utils ${LIBS} distributions)


#add_executable(pp utils/post_processor.cc parsers.cc model.cc)
#target_link_libraries(pp ${Boost_LIBRARIES} ${Bamtools_LIBRARIES})
#link_directories("${Bamtools_LIBRARY_DIR}/bamtools/")
#target_link_libraries(pp ${Boost_LIBRARIES} bamtools bamtools-utils )


#add_executable(freq_seq utils/freq_seq.cc parsers.cc)
#target_link_libraries(pp bamtools bamtools-utils)
