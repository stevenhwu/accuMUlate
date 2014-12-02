/*
 * MutationProb.h
 *
 *  Created on: Nov 8, 2014
 *      Author: Steven Wu
 */

#ifndef MUTATIONPROB_H_
#define MUTATIONPROB_H_

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED(func) func
#endif

#include "model.h"
#include "Lookup.h"



typedef std::array<double, 10> Array10D;
typedef std::array<double, 4> Array4D;

struct MutationRate{
    double prob;
    double one_minus_p;
};


class MutationProb {
public:
    static double CalculateBeta0(Array4D freq);
    static double CalculateExpBeta(double mu, double beta0);
    static double ConvertExpBetaToMu(double exp_beta, double beta0);


public:
	MutationProb(const ModelParams &model_params);
	~MutationProb();

    void UpdateMu(double mu);
    double ConvertExpBetaToMu(double exp_beta);

    MutationRate GetMutationRate();
    Array4D GetFrequencyPrior();
    Array10D GetAncestorPrior();
    double GetExpBeta();
    double GetMu0();
    double GetBeta0();


private:

    Array4D frequency_prior;
    Array10D ancestor_prior;
    MutationRate mutation_rate;

    double mu0;
    double beta0;
    double exp_beta; //exp_beta = exp(-beta0*mu0)

    double CalculateExpBeta();

};
#endif /* MUTATIONPROB_H_ */
