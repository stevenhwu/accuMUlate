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

class EmSummaryStatBinomial : public EmSummaryStat {

public:


    EmSummaryStatBinomial();

    void SetStats(double stat_same0, double stat_diff0);

    virtual void print();

    void Reset();

    void UpdateSumWithProportion(double &d, EmSummaryStatBinomial mutation);


    virtual ~EmSummaryStatBinomial() {
    }

//private:
    double stat_diff = 0;

};


#endif //EM_SUMMARY_STAT_BINOMIAL_H_
