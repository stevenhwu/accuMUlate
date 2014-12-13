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


    virtual void UpdateSummaryStat(double &prob, std::vector<double> &temp_stat);
    void UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat);

    virtual void UpdateEmModel(EmModel *em_model);



public:

    std::vector<int> data;
    int count;
    int count_negative;
    int total_count;

    double binomial_prob;

};


#endif //EM_DATA_BINARY_H_
