/*
 * pileup_utils.h
 *
 *  Created on: 12/13/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef PILEUP_UTILS_H_
#define PILEUP_UTILS_H_

#include <iostream>
#include <vector>
#include <boost/program_options.hpp>
#include <time.h>


//#include "api/BamReader.h"
#include "local_bamtools/bamtools_pileup_engine.h"

#include "local_bamtools/bamtools_fasta.h"
#include <mutations/model.h>
#include "VariantVisitor.h"
#include "variant_visitor_two.h"

#include "genome_data_stream.h"
#include "boost_input_utils.h"

//#include "parsers.h"
//#include "sequence_prob.h"

namespace PileupUtils {

    extern uint64_t print_count;
    extern float fix_counter;
    extern float time_scaler;



    void CreatePileupAlignment(boost::program_options::variables_map &vm, GenomeData &genome_data, int index);

    void CreatePileupV2(boost::program_options::variables_map &vm, GenomeData &genome_data, RefVector &references,
                        SamHeader &header, LocalBamToolsUtils::Fasta &reference_genome,
                        LocalBamToolsUtils::PileupEngine &pileup);

    void CreatePileupV1(boost::program_options::variables_map &vm, GenomeData &genome_data, RefVector &references,
                        SamHeader &header, LocalBamToolsUtils::Fasta &reference_genome,
                        LocalBamToolsUtils::PileupEngine &pileup);

    void WriteGenomeDataToBinary(std::string file_name, GenomeData &base_counts) ;

    void ReadGenomeDataFromBinary(std::string file_name, GenomeData &genome_data) ;

    void VerboseAlignmentInfo(clock_t &time_stored, uint64_t ali_counter);

    void SummariseReadsData(GenomeData base_counts);
};


#endif //PILEUP_UTILS_H_
