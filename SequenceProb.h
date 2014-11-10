/*
 * SequencingProb.h
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */

#ifndef SEQUENCINGPROB_H_
#define SEQUENCINGPROB_H_

#include <vector>
#include <model.h>
#include "MutationProb.h"
#include "distributions/DirichletMultinomialDistribution.h"


class SequenceProb{
	DiploidProbs anc_genotypes;
	DiploidProbs pop_genotypes;
	std::vector<HaploidProbs> allHap;
	ModelParams params;

	MutationMatrix non_mutation;
	MutationMatrix mutation;

	double likelihood;

	ModelInput data;

public:
	SequenceProb(){};
	SequenceProb(const ModelParams &params, const ModelInput site_data, MutationProb muProb);
	~SequenceProb();

//	void UpdateMu();
//	void UpdateMu(double mu);
	void UpdateLikelihood();
	void UpdateMuProb(MutationProb muProb);

	double GetLikelihood();

	ModelInput GetData();
protected:
	DiploidProbs DiploidPopulation(int ref_allele);
	HaploidProbs HaploidSequencing(ReadData data);
	DiploidProbs DiploidSequencing(ReadData data);

//	MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut);
//	MutationMatrix MutationAccumulation2(bool and_mut);



};

#endif /* SEQUENCINGPROB_H_ */
