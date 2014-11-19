#pragma once
#ifndef __JC69_H_
#define __JC69_H_

#include <array>
#include "EvolutionModel.h"

class JC69 : public EvolutionModel {


    const double FOUR_THIRD = 4.0 / 3.0;
    const double QUARTER = 1.0 / 4.0;
    const double THREE_QUARTER = 3.0 / 4.0;

public:
    JC69(double mu);

    ~JC69();

//    virtual double GetConditionalProb(int i);
    virtual std::array<double, 3> GetConditionalProbSpecial();

private:
    double exp_four_third_mu;
    std::array<double, 3> conditional_prob;

    virtual void UpdateConditionalProb();

    void computeExpFourThirdMu();
    double probNotEqual();
    double probOneEqual();
    double probThreeEqual();


};


#endif //__JC69_H_
