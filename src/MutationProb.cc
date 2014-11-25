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
    for (int i = 0; i < 4; ++i) {
        frequency_prior[i] = (i+1)/10.0;
        frequency_prior[i] = model_params.nuc_freq[i];
    }
    for (int i = 0; i < 4; ++i) {

        beta0 -= (frequency_prior[i] * frequency_prior[i]);
//        beta0 += params.nuc_freq[i] * params.nuc_freq[i];
        for (int j = i; j < 4; ++j) {
            int index10 = LookupTable::index_converter_16_to_10[i][j];
            ancestor_prior[index10] = frequency_prior[i] * frequency_prior[j];
            if(i != j){
                ancestor_prior[index10] *= 2; //Count both AC and CA
            }
        }
    }
    beta0 = 1.0 / beta0;

    UpdateMu(model_params.mutation_rate);
}

MutationProb::~MutationProb() {
}


void MutationProb::UpdateMu(double mu0) {
    this->mu0 = mu0;
    CalculateBeta();

	mutation_rate.mu = 1-exp_beta;
    mutation_rate.one_minus_mu = 1-mutation_rate.mu;


}

double MutationProb::CalculateBeta(){

    exp_beta = exp(-beta0 * mu0);
    return exp_beta;
}

//MutationMatrix MutationProb::GetMutation() {
//
//	return mutation;
//}
//
//MutationMatrix MutationProb::GetNonMutation() {
//	return non_mutation;
//}


MutationRate MutationProb::GetMutationRate() {
    return mutation_rate;
}

double MutationProb::GetBeta() {
    return exp_beta;
}

Array4D MutationProb::GetFrequencyPrior() {
    return frequency_prior;
}

Array10D MutationProb::GetAncestorPrior() {
    return ancestor_prior;
}

double MutationProb::GetMu0() {
    return mu0;
}

double MutationProb::GetBeta0() {
    return beta0;
}