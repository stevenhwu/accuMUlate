/*
 * em_model_mutation.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_MODEL_MUTATION_V1_H_
#define EM_MODEL_MUTATION_V1_H_


#include "evolution_models/EvolutionModel.h"
#include "em_model.h"

class [[deprecated]] EmModelMutationV1 : public EmModel {

public:


    EmModelMutationV1(EvolutionModel &evo_model0);

    EmModelMutationV1(const EmModelMutationV1 &em_model);

    virtual ~EmModelMutationV1() {}

    virtual void UpdateParameter(double param);

    virtual void GetParameterInfo();

    virtual void UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat, double &log_likelihood_scaler);

    virtual size_t GetDataCount();

    EvolutionModel * GetEvoModel();



protected:
    EvolutionModel *evo_model;



};


#endif //EM_MODEL_MUTATION_V1_H_
