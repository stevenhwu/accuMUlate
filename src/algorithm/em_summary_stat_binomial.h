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
    static const int EM_SUMMARY_STAT_BINOMIAL_STATS_COUNT = 2;
    /*
    stat[0] = stat_same;
    stat[1] = stat_diff;
     */

    EmSummaryStatBinomial();

    virtual ~EmSummaryStatBinomial(){}

    virtual double MaximiseStats();
};


#endif //EM_SUMMARY_STAT_BINOMIAL_H_
