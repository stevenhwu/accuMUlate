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


    virtual void UpdateParameter(double param) = 0; //TODO: implement multi parameters later

    virtual void GetParameterInfo() = 0;

    virtual size_t GetDataCount() = 0;


    virtual void UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat) =0;
};


#endif //EM_MODEL_H_
