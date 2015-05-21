add_library(LibIoData
    local_bamtools/bamtools_fasta.cpp
    local_bamtools/bamtools_pileup_engine.cpp
    boost_input_utils.cc
    genome_data_stream.cc
    parsers.cc
    pileup_utils.cc
    variant_visitor_two.cc
    VariantVisitor.cc
)
#target_link_libraries(LibIoData ${LIBS} )
target_link_libraries(LibIoData ${Bamtools_LIBRARIES})

