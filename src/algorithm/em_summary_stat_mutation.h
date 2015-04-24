/*
 * em_summary_stat_mutation.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_SUMMARY_STAT_MUTATION_V1_H_
#define EM_SUMMARY_STAT_MUTATION_V1_H_


#include "em_summary_stat_binomial.h"


class EmSummaryStatMutation : public EmSummaryStatBinomial{


public:
    static const int EM_SUMMARY_STAT_MUTATION_STATS_COUNT = 2;
    EmSummaryStatMutation() : EmSummaryStatBinomial(){}

    virtual ~EmSummaryStatMutation() {}


    virtual double MaximiseStats() override;
};


#endif //EM_SUMMARY_STAT_MUTATION_V1_H_


