#pragma once
#ifndef __F81_H_
#define __F81_H_


#include <model.h>
#include "EvolutionModel.h"

class F81 : public EvolutionModel{
//    double mu;

public:
    F81(double mu, vector<double> params);
    ~F81();

protected:
    double beta0;
    vector<double> freqs;

    virtual void UpdateTransitionMatrix();



private:
    double exp_beta;

};


#endif //__F81_H_
