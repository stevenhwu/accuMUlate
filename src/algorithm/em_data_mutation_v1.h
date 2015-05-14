/*
 * em_data_mutation.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_DATA_MUTATION_V1_H
#define EM_DATA_MUTATION_V1_H_


#include <mutations/site_genotypes.h>
#include <mutations/site_prob.h>
#include "em_data.h"
#include "em_model.h"
#include "em_model_mutation_v1.h"




class EmDataMutationV1 : public EmData {
protected:
//    SiteProb site;

public:
    EmDataMutationV1(SiteGenotypes &sequence_prob, EvolutionModel &evo_model);

    SiteProb site;

    EmDataMutationV1() {}
    EmDataMutationV1(SequenceProb &sequence_prob, EvolutionModel &evo_model);
    virtual  ~EmDataMutationV1();

    virtual void UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat);

    virtual void UpdateSummaryStat(double &prob, std::vector<double> &temp_stat);

    virtual void UpdateEmModel(EmModel *em_model);
};


#endif //EM_DATA_MUTATION_V1_H
