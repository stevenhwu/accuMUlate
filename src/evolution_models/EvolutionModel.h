#include <model.h>

#pragma once
#ifndef __EvolutionModel_H_
#define __EvolutionModel_H_


class EvolutionModel {


public:

//    EvolutionModel(double mu);

    EvolutionModel(MutationProb mu_prob);

    virtual ~EvolutionModel();

    void UpdateMu(double mu0);

    MutationMatrix GetTranstionMatirxAToD();

    MutationRate GetMutationRate();

    MutationProb GetMutationProb();

protected:
    double mu;
    MutationProb mu_prob;
    MutationMatrix transition_matrix_a_to_d;

    virtual void UpdateTransitionMatrix() = 0;


};


#endif //__EvolutionModel_H_
