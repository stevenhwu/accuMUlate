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

    EmModelMutation(const EmModelMutation &em_model);

    virtual ~EmModelMutation() {
    }

    virtual void UpdateParameter(double param);

    void GetParameterInfo();

    EvolutionModel * GetEvoModel();
//    EvolutionModel * GetEvoModel() const;


protected:
    EvolutionModel *evo_model;

//
//public:
//    virtual EmModelMutation *GetModel();
//
//    std::unique_ptr<EvolutionModel> CopyEvoModel() const;
//    EvolutionModel* CopyEvoModel2()const ;
};


#endif //EM_MODEL_MUTATION_H_
