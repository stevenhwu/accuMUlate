/*
 * em_model_mutation.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_MODEL_MUTATION_H_
#define EM_MODEL_MUTATION_H_


#include "evolution_models/EvolutionModel.h"
#include "em_model.h"

class EmModelMutation : public EmModel {

public:
    EmModelMutation(EvolutionModel &evo_model0);

    EmModelMutation(EmModelMutation &em_model);

    virtual ~EmModelMutation() {
    }

    virtual void UpdateParameter(double param);

    void GetParameterInfo();

    EvolutionModel *GetEvoModel() ;

protected:
    EvolutionModel *evo_model;




};


#endif //EM_MODEL_MUTATION_H_
