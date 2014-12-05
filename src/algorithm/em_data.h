/*
 * em_data.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_DADA_H_
#define EM_DADA_H_


#include <iostream>
#include "em_summary_stat.h"
#include "em_model.h"


class EmData {

public:

    virtual ~EmData() {
    };

    virtual void UpdateSummaryStat(double prob, EmSummaryStat &summaryStat) = 0;
//    virtual void UpdateEmModel(EmModel &em_model) = 0;
    virtual void UpdateEmModel(EmModel *em_model) = 0;

    virtual void Test(double num) = 0;
};


#endif //EM_DADA_H_
