#pragma once
#ifndef __F81_H_
#define __F81_H_



#include "EvolutionModel.h"

class F81 : public EvolutionModel{


public:
//    F81(double mu, Array4D freqs);
//    F81(double mu, vector<double> freqs); //Legacy
    F81(MutationProb &prob);

    F81(double mu);

    ~F81();
    F81(const F81&); // copy ctor


    virtual std::unique_ptr<EvolutionModel> Clone() const;
    virtual EvolutionModel *Clone2() const;
    virtual void UpdateOneMinusExpBeta(double exp_beta);

    double m1();
    double m2();
protected:


    virtual void UpdateTransitionMatrix();



private:
    double exp_beta;
    double beta0;
    Array4D freqs;

    Eigen::Matrix4d m;


    void OtherTransitionMatrix();
};


#endif //__F81_H_
