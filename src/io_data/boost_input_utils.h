/*
 * boost_input_utils.h
 *
 *  Created on: 12/12/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef BOOST_INPUT_UTILS_H_
#define BOOST_INPUT_UTILS_H_


#include <boost/program_options/variables_map.hpp>

#include <boost/program_options.hpp>
#include <time.h>


#include "io_data/parsers.h"
//#include "api/BamReader.h"
//#include "utils/bamtools_pileup_engine.h"
//#include "utils/bamtools_fasta.h"

//#include "third-party/bamtools/src/api/BamReader.h"
//#include "third-party/bamtools/src/utils/bamtools_pileup_engine.h"
//#include "third-party/bamtools/src/utils/bamtools_fasta.h"

#include "api/BamReader.h"
//#include "utils/bamtools_fasta.h"
//#include <api/BamReader.h>
//#include <utils/bamtools_fasta.h>

#include "local_bamtools/bamtools_fasta.h"
#include <mutations/model.h>
#include <sys/stat.h>



namespace BoostUtils {



//    namespace po = boost::program_options;

    void ParseCommandLinkeInput(int argc, char **argv, boost::program_options::variables_map &vm);

    void ExtractInputVariables(boost::program_options::variables_map &vm, GenomeData &genome_data,
            BamTools::BamReader &experiment, BamTools::RefVector &references,
            BamTools::SamHeader &header, LocalBamToolsUtils::Fasta &reference_genome);

    ModelParams CreateModelParams(boost::program_options::variables_map variables_map);


}

#endif //BOOST_INPUT_UTILS_H_
