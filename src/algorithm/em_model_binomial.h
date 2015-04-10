/*
 * em_model_binary.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_MODEL_BINARY_H_
#define EM_MODEL_BINARY_H_

#include <stddef.h>
#include "em_model.h"


class EmModelBinomial : public EmModel {

public:

    EmModelBinomial(int n0, double prob0);

    virtual ~EmModelBinomial(){}

    virtual void UpdateParameter(double param);

    virtual size_t GetDataCount();

    virtual void UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat, double &log_likelihood_scaler);

//    virtual void GetParameterInfo();

    double GetParameter();

protected:
    int n;
    double prob;

};


#endif //EM_MODEL_BINARY_H_
