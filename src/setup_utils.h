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
#include <algorithm/em_algorithm_thread_mutation.h>

GenomeData GetGenomeData(boost::program_options::variables_map variables_map);

void SummariseRealData(GenomeData &genome_data );

void SummariseRealData(boost::program_options::variables_map map);

void CustomFilterGenomeData(GenomeData &genome_data);


void RemoveDescendantRDV(ReadDataVector &read_vector, int d);

void RemoveSiteGenomeData(GenomeData &genome_data, int g);


class setup_utils {

};


#endif //ACCUMULATE_SETUP_UTILS_H
