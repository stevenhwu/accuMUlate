

SET(SOURCE_FILES ${SOURCE_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/em_algorithm.cc
    #${CMAKE_CURRENT_SOURCE_DIR}/em_data_old.cc
    #${CMAKE_CURRENT_SOURCE_DIR}/em
    ${CMAKE_CURRENT_SOURCE_DIR}/em_data_mutation.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_model_mutation.cc

    ${CMAKE_CURRENT_SOURCE_DIR}/em_summary_stat_binomial.cc
#    ${CMAKE_CURRENT_SOURCE_DIR}/ em_summary_stat_binomial.cc
    #PARENT_SCOPE
    CACHE INTERNAL "Global"
    )
#list (APPEND SOURCE_LIST ${SOURCE_FILES})
message("Algorithm: ${SOURCE_FILES} ")
