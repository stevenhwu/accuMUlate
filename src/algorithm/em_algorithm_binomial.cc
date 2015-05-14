/*
 * em_algorithm_binomial.cc
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */


#include <stddef.h>
#include "em_algorithm_binomial.h"
#include "em_summary_stat_binomial.h"

EmAlgorithmBinomial::EmAlgorithmBinomial(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr0, EmModelBinomial &em_model0)
        : EmAlgorithm(num_category0, data_ptr0, em_model0) {

   for (size_t i = 0; i < num_category; ++i) {
        em_model.emplace_back(new EmModelBinomial(em_model0));
    }

    InitWithData();
}

void EmAlgorithmBinomial::InitialiseSummaryStat() {

//    em_stat_local_single = std::unique_ptr<EmSummaryStat>(new EmSummaryStatBinomial());
//    em_stat_local_single->Print();
    temp_stats = std::vector<std::vector<double>>(num_category);
    for (size_t i = 0; i < num_category; ++i) {
        all_em_stats.emplace_back(new EmSummaryStatBinomial());
        temp_stats[i] = std::vector<double>(all_em_stats[i]->GetStatCount());
    }

}

void EmAlgorithmBinomial::RunEM() {

    size_t i = 0;
    bool isConverged = true;
    while(isConverged){

        ExpectationStepModel();
        MaximizationStep();

        isConverged = EmStoppingCriteria(i);
        i++;
    }


}



void EmAlgorithmBinomial::InitialiseParameters() {
    double lower_bound = 0.5;
    double upper_bound = 0.6;

    if (num_category == 2) {
        parameters = {upper_bound, lower_bound};
    }
    else {
        std::cout << "Not yet implemented for more than 2 categories" << std::endl;
        exit(222);
        //TODO: Should throw exception instead of exit, this will do for now
    }

}

void EmAlgorithmBinomial::ExpectationStepCustom(size_t data_index, size_t category_index, double &sum_prob, std::vector<double> &temp_stat) {


    em_data_ptr->at(data_index)->UpdateEmModel( em_model[category_index].get() );
    em_data_ptr->at(data_index)->UpdateSummaryStat(sum_prob, temp_stat);


}
