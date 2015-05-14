/*
 * mutation_prob.h
 *
 *  Created on: Nov 8, 2014
 *      Author: Steven Wu
 */
#pragma once
#ifndef MUTATION_PROB_H_
#define MUTATION_PROB_H_

//#ifdef __GNUC__
//#define DEPRECATED(func) func __attribute__ ((deprecated))
//#elif defined(_MSC_VER)
//#define DEPRECATED(func) __declspec(deprecated) func
//#else
//#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
//#define DEPRECATED(func) func
//#endif

#include <mutations/model.h>
#include "lookup.h"



typedef std::array<double, 10> Array10D;
typedef std::array<double, 4> Array4D;


class MutationProb {
public:
    static double CalculateBeta0(Array4D &freq);
    static double CalculateExpBeta(double mu, double beta0);
    static double ConvertOneMinusExpBetaToMu(double one_minus_exp_beta, double beta0);


public:
    MutationProb();
    explicit MutationProb(double mu);
//    MutationProb(const MutationProb& self);
    MutationProb(const ModelParams &model_params);
	~MutationProb();

    void UpdateMu(double mu);
    double ConvertOneMinusExpBetaToMu(double one_minus_exp_beta);

    double GetMutationRate();
    Array4D GetFrequencyPrior();
    Array10D GetAncestorPrior();
    double GetExpBeta();
    double GetMu();
    double GetBeta0();


private:
    Array4D frequency_prior;
    Array10D ancestor_prior;

    double mutation_rate;
    double mu;
    double beta0;
    double exp_beta; //exp_beta = exp(-beta0*mu0)

    double CalculateExpBeta();
    void CalculateAncestorPrior();

    static double ConvertExpBetaToMu(double exp_beta0, double beta0);
};
#endif /* MUTATION_PROB_H_ */
