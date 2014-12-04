/*
 * em_data_mutation.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */

#pragma once

#include "em_data.h"
#include "sequence_prob.h"
#include "site_prob.h"

#ifndef EM_DATA_MUTATION_H
#define EM_DATA_MUTATION_H_




class EmDataMutation : public EmData {

public:


    EmDataMutation() {};


    virtual  ~EmDataMutation();
    EmDataMutation(SequenceProb &sequence_prob, EvolutionModel &evo_model);

    virtual void UpdateLikelihood(double prob, EmSummaryStat &summaryStat);

    SiteProb site;
};


#endif //EM_DATA_MUTATION_H
