

SET(SOURCE_FILES ${SOURCE_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/em_algorithm.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_algorithm_thread.cc

    ${CMAKE_CURRENT_SOURCE_DIR}/em_algorithm_binomial.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_algorithm_mutation_v1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_algorithm_mutation.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_algorithm_thread_mutation.cc

    ${CMAKE_CURRENT_SOURCE_DIR}/em_data_binomial.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_data_mutation_v1.cc

    ${CMAKE_CURRENT_SOURCE_DIR}/em_model_binomial.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_model_mutation_v1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_model_mutation.cc

    ${CMAKE_CURRENT_SOURCE_DIR}/em_logger.cpp #algorithm/em_logger.h

    ${CMAKE_CURRENT_SOURCE_DIR}/em_summary_stat.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_summary_stat_binomial.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/em_summary_stat_mutation.cc

    #PARENT_SCOPE
    CACHE INTERNAL "Global"
    )
#list (APPEND SOURCE_LIST ${SOURCE_FILES})
#message("Algorithm: ${SOURCE_FILES} ")
