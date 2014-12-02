#pragma once
#ifndef __F81_H_
#define __F81_H_


#include <model.h>
#include <MutationProb.h>
#include "EvolutionModel.h"

class F81 : public EvolutionModel{
//    double mu;

public:
    F81(double mu, Array4D freqs);
    F81(double mu, vector<double> freqs); //Legacy
    ~F81();

protected:
    double beta0;
    Array4D freqs;

    virtual void UpdateTransitionMatrix();



private:
    double exp_beta;

};


#endif //__F81_H_
