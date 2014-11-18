/*
 * MutationProb.h
 *
 *  Created on: Nov 8, 2014
 *      Author: Steven Wu
 */

#ifndef MUTATIONPROB_H_
#define MUTATIONPROB_H_

#include "model.h"

class MutationProb {



	MutationMatrix non_mutation;
	MutationMatrix mutation;

	ModelParams params;
	int beta0;

public:
	MutationProb(const ModelParams &model_params);
	~MutationProb();

	void UpdateMu();
	void UpdateMu(double mu);
	MutationMatrix GetMutation();
	MutationMatrix GetNonMutation();

    double GetMu();

protected:

	MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut);

	MutationMatrix MutationAccumulation2(bool and_mut);



private:
	double CalculateBeta();

};
#endif /* MUTATIONPROB_H_ */
