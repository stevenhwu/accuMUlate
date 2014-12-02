#include <model.h>

#pragma once
#ifndef __EvolutionModel_H_
#define __EvolutionModel_H_


class EvolutionModel {


public:

    EvolutionModel(double mu);

    virtual ~EvolutionModel();

    void UpdateMu(double mu);

    MutationMatrix GetTranstionMatirxAToD();

protected:
    double mu;
    MutationMatrix transition_matrix_a_to_d;

    virtual void UpdateTransitionMatrix() = 0;


};


#endif //__EvolutionModel_H_
