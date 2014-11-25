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
//#include "SequenceProb.h"


typedef std::array<double, 10> Array10D;
typedef std::array<double, 4> Array4D;

struct MutationRate{
    double mu;
    double one_minus_mu;
};


class MutationProb {

//	MutationMatrix non_mutation;
//	MutationMatrix mutation;
//    ModelParams params;
    double mu0;
	double beta0;
    double exp_beta; //exp_beta = exp(-beta*mu)

public:
	MutationProb(const ModelParams &model_params);
	~MutationProb();

//	void UpdateMu();
	void UpdateMu(double mu);
//	MutationMatrix GetMutation();
//	MutationMatrix GetNonMutation();

    MutationRate GetMutationRate();
    double GetBeta();

    Array4D GetFrequencyPrior();
    Array10D GetAncestorPrior();

    double GetMu0();
    double GetBeta0();

protected:

    MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut) ;
    MutationMatrix MutationAccumulation2(bool and_mut) ;



private:
	double CalculateBeta();

    Array4D frequency_prior;
    Array10D ancestor_prior;
    MutationRate mutation_rate;


};
#endif /* MUTATIONPROB_H_ */
