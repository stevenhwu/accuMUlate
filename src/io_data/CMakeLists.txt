add_library(LibIoData
    bamtools_fasta.cpp
    boost_input_utils.cc
    genome_data_stream.cc
    parsers.cc
    pileup_utils.cc
    variant_visitor_two.cc
    VariantVisitor.cc
)
#target_link_libraries(LibIoData ${LIBS} )
target_link_libraries(LibIoData ${Bamtools_LIBRARIES})

