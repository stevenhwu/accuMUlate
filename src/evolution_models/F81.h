#pragma once
#ifndef __F81_H_
#define __F81_H_


#include <model.h>
#include "mutation_prob.h"
#include "EvolutionModel.h"

class F81 : public EvolutionModel{


public:
//    F81(double mu, Array4D freqs);
//    F81(double mu, vector<double> freqs); //Legacy
    F81(MutationProb prob);

    F81(double mu);

    ~F81();

protected:


    virtual void UpdateTransitionMatrix();



private:
    double exp_beta;
    double beta0;
    Array4D freqs;
};


#endif //__F81_H_
