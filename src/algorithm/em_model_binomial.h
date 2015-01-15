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

    virtual ~EmModelBinomial(){};

    void UpdateParameter(double param);

    double GetParameter();

    virtual void GetParameterInfo();

    int n;
    double prob;

    virtual size_t GetDataCount();
};


#endif //EM_MODEL_BINARY_H_
