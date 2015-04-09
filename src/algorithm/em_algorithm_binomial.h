/*
 * em_algorithm_binomial.h
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_ALGORITHM_BINOMIAL_H_
#define EM_ALGORITHM_BINOMIAL_H_

#include <stddef.h>
#include "em_algorithm.h"
class EmAlgorithmBinomial : public EmAlgorithm {

public:
    EmAlgorithmBinomial(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr0, EmModelBinomial &em_model0);

    EmAlgorithmBinomial(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr0, EmModel &em_model0);

    virtual void RunEM();


protected:

    virtual void ExpectationStepCustom(size_t data_index, size_t category_index, double &sum_prob, std::vector<double> &temp_stat);

    virtual void InitialiseSummaryStat();


    void InitialiseParameters();
};


#endif //EM_ALGORITHM_BINOMIAL_H_
