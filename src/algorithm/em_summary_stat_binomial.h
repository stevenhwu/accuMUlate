/*
 * em_summary_stat_binomial.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_SUMMARY_STAT_BINOMIAL_H_
#define EM_SUMMARY_STAT_BINOMIAL_H_


#include "em_summary_stat.h"
#include <vector>
class EmSummaryStatBinomial : public EmSummaryStat {

public:


    EmSummaryStatBinomial();

    virtual ~EmSummaryStatBinomial() {
    }


    void biTest(double a);
};


#endif //EM_SUMMARY_STAT_BINOMIAL_H_
