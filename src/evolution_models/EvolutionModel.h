#include <memory>
#include <mutations/mutation_prob.h>


#include <mutations/model.h>

#pragma once
#ifndef __EvolutionModel_H_
#define __EvolutionModel_H_


class EvolutionModel {


public:

//    EvolutionModel(double mu);
//    EvolutionModel(const EvolutionModel & evolution_model);

    EvolutionModel(MutationProb mu_prob);



    virtual ~EvolutionModel();

    EvolutionModel() {
    }

    void UpdateMu(double mu0);

    MutationMatrix GetTranstionMatirxAToD();

    double GetMutationRate();

    MutationProb GetMutationProb();
    MutationProb GetMutationProb() const;

    virtual void UpdateOneMinusExpBeta(double d) = 0;//FIXME, clarify expBeta or 1-expBeta, should be 1-expBeta here

    virtual std::unique_ptr<EvolutionModel> Clone() const = 0;
    virtual EvolutionModel *Clone2() const = 0;

protected:
    double mu;
    MutationProb mu_prob;
    MutationMatrix transition_matrix_a_to_d;

    virtual void UpdateTransitionMatrix() = 0;

};


#endif //__EvolutionModel_H_
