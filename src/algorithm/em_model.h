/*
 * em_model.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_MODEL_H_
#define EM_MODEL_H_

#include <vector>
class EmModel {
public:


    virtual ~EmModel() {
    }

    virtual void UpdateParameter(double param) = 0; //TODO: implement multi parameters later

    virtual void UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat, double &log_likelihood_scaler) =0;

    virtual size_t GetDataCount() = 0;

    [[deprecated]] virtual void GetParameterInfo(){}

};


#endif //EM_MODEL_H_
