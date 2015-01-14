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

class EmModelMutationV1 : public EmModel {

public:
    EmModelMutationV1(EvolutionModel &evo_model0);

    EmModelMutationV1(const EmModelMutationV1 &em_model);

    virtual ~EmModelMutationV1() {
    }

    virtual void UpdateParameter(double param);

    void GetParameterInfo();

    EvolutionModel * GetEvoModel();
//    EvolutionModel * GetEvoModel() const;


protected:
    EvolutionModel *evo_model;

//
//public:
//    virtual EmModelMutationV1 *GetModel();
//
//    std::unique_ptr<EvolutionModel> CopyEvoModel() const;
//    EvolutionModel* CopyEvoModel2()const ;
};


#endif //EM_MODEL_MUTATION_V1_H_
