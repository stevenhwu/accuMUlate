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
#include "em_model_mutation.h"
#include "em_model_binomial.h"


#ifndef EM_DATA_MUTATION_H
#define EM_DATA_MUTATION_H_




class EmDataMutation : public EmData {
protected:
//    SiteProb site;

public:
    SiteProb site;

    EmDataMutation() {};
    EmDataMutation(SequenceProb &sequence_prob, EvolutionModel &evo_model);
    virtual  ~EmDataMutation();

    virtual void UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat);

//    virtual void UpdateEmModel(EmModel &em_model);





    virtual void UpdateEmModel(EmModel *em_model);


    virtual void UpdateSummaryStat(double &prob, std::vector<double> &temp_stat);
};


#endif //EM_DATA_MUTATION_H
