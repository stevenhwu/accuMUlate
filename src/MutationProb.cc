/*
 * MutationProb.cpp
 *
 *  Created on: Nov 8, 2014
 *      Author: Steven Wu
 */

#include "MutationProb.h"
#include <cmath>
#include <iostream>

#include "SequenceProb.h"

MutationProb::MutationProb(const ModelParams &model_params) {

	beta0 = 1.0;
	for (auto d : model_params.nuc_freq) {
		beta0 -= d * d;
	}
	beta0 = 1.0 / beta0;

    for (int i = 0; i < 4; ++i) {
        frequency_prior[i] = model_params.nuc_freq[i];
//        beta0 += params.nuc_freq[i] * params.nuc_freq[i];
        for (int j = i; j < 4; ++j) {
            int index10 = LookupTable::index_converter_16_to_10[i][j];
            ancestor_prior[index10] = model_params.nuc_freq[i] * model_params.nuc_freq[j];
            if(i != j){
                ancestor_prior[index10] *= 2; //Count both AC and CA
            }
        }
    }

    UpdateMu(model_params.mutation_rate);
}

MutationProb::~MutationProb() {
}


void MutationProb::UpdateMu(double mu0) {
	mutation_rate.mu = mu0;
    mutation_rate.one_minus_mu = 1 - mutation_rate.mu;
    CalculateBeta();

}

MutationMatrix MutationProb::GetMutation() {

	return mutation;
}

MutationMatrix MutationProb::GetNonMutation() {
	return non_mutation;
}

double MutationProb::CalculateBeta(){

	beta = exp(-beta0 * mutation_rate.mu);
	return beta;
}

MutationRate MutationProb::GetMutationRate() {
    return mutation_rate;
}

double MutationProb::GetBeta() {
    return beta;
}

Array4D MutationProb::GetFrequencyPrior() {
    return frequency_prior;
}

Array10D MutationProb::GetAncestorPrior() {
    return ancestor_prior;
}
