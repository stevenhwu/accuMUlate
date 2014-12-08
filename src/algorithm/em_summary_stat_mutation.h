/*
 * em_summary_stat_mutation.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_SUMMARY_STAT_MUTATION_H_
#define EM_SUMMARY_STAT_MUTATION_H_

#include "em_summary_stat_binomial.h"
#include "em_summary_stat.h"

class EmSummaryStatMutation : public EmSummaryStatBinomial{

public:
    EmSummaryStatMutation() {
    }

    virtual ~EmSummaryStatMutation() {
    }


    virtual double MaximiseStats() override;
};


#endif //EM_SUMMARY_STAT_MUTATION_H_


