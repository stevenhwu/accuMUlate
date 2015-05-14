/*
 * em_model_mutation.h
 *
 *  Created on: 1/14/15
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_MODEL_MUTATION_H_
#define EM_MODEL_MUTATION_H_


#include <mutations/model.h>
#include <mutations/mutation_model.h>
#include "em_model.h"

class EmModelMutation : public EmModel {

public:
    EmModelMutation(MutationModel &mutation_model);

//    EmModelMutation(const MutationModel &mutation_model);

    virtual ~EmModelMutation() {
    }

    virtual void UpdateParameter(double param);

    virtual void UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat, double &log_likelihood_scaler);

    virtual size_t GetDataCount();

    virtual void GetParameterInfo();

//    EvolutionModel * GetEvoModel();
//    EvolutionModel * GetEvoModel() const;


protected:
    MutationModel mutation_model;


};


#endif //EM_MODEL_MUTATION_H_
