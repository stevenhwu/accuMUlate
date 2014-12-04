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

class EmData {

public:

    virtual ~EmData() {
    };

    virtual void UpdateLikelihood(double prob, EmSummaryStat &summaryStat) = 0;
};


#endif //EM_DADA_H_
