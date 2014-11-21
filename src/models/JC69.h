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

    std::array<double, 3> GetConditionalProbSpecial();


protected:
    virtual void UpdateTransitionMatrix();

private:
    double exp_four_third_mu;
    std::array<double, 3> conditional_prob_special;

    virtual void UpdateConditionalProbSpecial();

    void computeExpFourThirdMu();
    double probNotEqual();
    double probOneEqual();
    double probThreeEqual();


};


#endif //__JC69_H_
