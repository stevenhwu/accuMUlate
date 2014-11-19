#include <model.h>

#pragma once
#ifndef __EvolutionModel_H_
#define __EvolutionModel_H_



class EvolutionModel {


public:

    EvolutionModel(double mu);

    virtual ~EvolutionModel();

    void UpdateMu(double mu);
    virtual MutationMatrix GetConditionalProb() = 0;

protected:
    double mu;

    virtual void UpdateConditionalProb() = 0;

};


#endif //__EvolutionModel_H_
