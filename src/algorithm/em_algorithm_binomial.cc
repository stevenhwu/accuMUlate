/*
 * em_algorithm_binomial.cc
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */


#include "em_algorithm_binomial.h"
#include "em_summary_stat_binomial.h"

EmAlgorithmBinomial::EmAlgorithmBinomial(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr0, EmModelBinomial &em_model0)
        : EmAlgorithm(num_category0, data_ptr0, em_model0) {

   for (int i = 0; i < num_category; ++i) {
        em_model.emplace_back(new EmModelBinomial(em_model0));
    }

    Init();
    em_count = 100;
}

void EmAlgorithmBinomial::InitialiseSummaryStat() {

    em_stat_local_single = std::unique_ptr<EmSummaryStat>(new EmSummaryStatBinomial());
    em_stat_local_single->print();
    temp_stats = std::vector<std::vector<double>>(num_category);
    for (size_t i = 0; i < num_category; ++i) {
        all_em_stats.emplace_back(new EmSummaryStatBinomial());
        temp_stats[i] = vector<double>(em_stat_local_single->GetStatCount());
    }

}

void EmAlgorithmBinomial::Run() {

//    em_stat_local_single->print();
    for (size_t i = 0; i < em_count; ++i) {
        ExpectationStep();
        MaximizationStep();
    }

}



void EmAlgorithmBinomial::InitialiseParameters() {
    double lower_bound = 0.1;
    double upper_bound = 0.9;
    parameters = std::vector<double>(num_category);
    if (num_category == 2) {
        parameters = {upper_bound, lower_bound};
    }
    else {
        std::cout << "Not yet implemented for more than 2 categories" << std::endl;
        exit(222);
        //TODO: Should throw exception instead of exit, this will do for now
    }

}

