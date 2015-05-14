

SET(SOURCE_FILES ${SOURCE_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/EvolutionModel.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/JC69.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/F81.cc
    #PARENT_SCOPE
    CACHE INTERNAL "Global"
    )
add_library(LibEvolutionModels EvolutionModel.cc F81.cc JC69.cc)
#list (APPEND SOURCE_LIST ${SOURCE_FILES})
#message("MODEL: ${SOURCE_FILES} ")
