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
#include "em_model.h"

#ifndef EM_DATA_MUTATION_H
#define EM_DATA_MUTATION_H_




class EmDataMutation : public EmData {

public:
    EmDataMutation() {};
    EmDataMutation(SequenceProb &sequence_prob, EvolutionModel &evo_model);
    virtual  ~EmDataMutation();

    virtual void UpdateSummaryStat(double prob, EmSummaryStat &summaryStat);

//    virtual void UpdateEmModel(EmModel &em_model);




    SiteProb site;

    virtual void Test(double num);

    void UpdateEmModel(EmModel *em_model);
};


#endif //EM_DATA_MUTATION_H
