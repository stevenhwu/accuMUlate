/*
 * em_algorithm_binomial.h
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_ALGORITHM_BINOMIAL_H_
#define EM_ALGORITHM_BINOMIAL_H_

#include "em_algorithm.h"
class EmAlgorithmBinomial : public EmAlgorithm {

public:
    EmAlgorithmBinomial(int num_category0, vector<unique_ptr<EmData>> &data_ptr0, EmModelBinomial &em_model0);

    virtual void Run();
       EmAlgorithmBinomial(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr0, EmModel &em_model0);

protected:
    virtual void InitialiseSummaryStat();


    void InitialiseParameters();
};


#endif //EM_ALGORITHM_BINOMIAL_H_