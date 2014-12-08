/*
 * em_data_binary.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_DATA_BINARY_H_
#define EM_DATA_BINARY_H_

#include "em_data.h"
#include "em_model_binomial.h"
#include <vector>

class EmDataBinomial : public EmData {

public:
    EmDataBinomial(std::vector<int> input);



    EmDataBinomial(int total, int count0);

    ~EmDataBinomial();


    std::vector<int> data;

    int count;

    virtual void UpdateEmModel(std::unique_ptr<EmModelBinomial> &em_model);

    virtual void UpdateEmModel(std::unique_ptr<EmModelMutation> &em_model);

    void UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat);

    void UpdateEmModel(EmModel *em_model);

    unsigned long total_count;
    double binomial_prob;
    unsigned long count_negative;
};


#endif //EM_DATA_BINARY_H_
