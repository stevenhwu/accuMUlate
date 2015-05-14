//
// Created by steven on 5/13/15.
//

#ifndef ACCUMULATE_SETUP_UTILS_H
#define ACCUMULATE_SETUP_UTILS_H


#include <iostream>
//#include <map>
//#include <vector>
#include <boost/program_options.hpp>
//#include <algorithm/em_algorithm_thread_mutation.h>
//#include <algorithm/em_algorithm_mutation.h>
#include <io_data/boost_input_utils.h>
#include <io_data/pileup_utils.h>

void SummariseReadsData(GenomeData base_counts);

void SummariseRealData(GenomeData &genome_data );


GenomeData getGenomeData(boost::program_options::variables_map variables_map);

void SummariseRealData(boost::program_options::variables_map map);

void PostFilterGenomeData(GenomeData &genome_data);


class setup_utils {

};


#endif //ACCUMULATE_SETUP_UTILS_H
