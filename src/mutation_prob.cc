/*
 * mutation_prob.cpp
 *
 *  Created on: Nov 8, 2014
 *      Author: Steven Wu
 */

#include "mutation_prob.h"
#include <cmath>
#include <iostream>

MutationProb::MutationProb(const ModelParams &model_params) {

    for (int i = 0; i < 4; ++i) {
        frequency_prior[i] = model_params.nuc_freq[i];
    }

    CalculateAncestorPrior();
    beta0 = CalculateBeta0(frequency_prior);
    UpdateMu(model_params.mutation_rate);
}


MutationProb::MutationProb(double mu) {

    for (int i = 0; i < 4; ++i) {
        frequency_prior[i] = 0.25;
    }
    CalculateAncestorPrior();
    beta0 = CalculateBeta0(frequency_prior);
    UpdateMu(mu);

}



MutationProb::~MutationProb() {
}


void MutationProb::CalculateAncestorPrior() {
    for (int i = 0; i < 4; ++i) {
        for (int j = i; j < 4; ++j) {
            int index10 = LookupTable::index_converter_4_4_to_10[i][j];
            ancestor_prior[index10] = frequency_prior[i] * frequency_prior[j];
            if(i != j){
                ancestor_prior[index10] *= 2; //Count both AC and CA
            }
        }
    }
}

void MutationProb::UpdateMu(double mu0) {
    this->mu0 = mu0;
    CalculateExpBeta();

	mutation_rate.prob = 1-exp_beta;
    mutation_rate.one_minus_p = 1-mutation_rate.prob;

//    cout << "MuP:\t" << this->mu0 << "\t" << mutation_rate.prob << endl;
}



double MutationProb::CalculateExpBeta(){
    exp_beta = CalculateExpBeta(mu0, beta0);
    return exp_beta;
}

double MutationProb::ConvertExpBetaToMu(double exp_beta) {
    return ConvertExpBetaToMu(exp_beta, beta0);
}

MutationRate MutationProb::GetMutationRate() {
    return mutation_rate;
}

double MutationProb::GetExpBeta() {
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


//////////////////////////////////////////////////
//// Static methods
//////////////////////////////////////////////////
double MutationProb::CalculateExpBeta(double mu, double beta0){
    double exp_beta = exp(-beta0 * mu);
    return exp_beta;
}


double MutationProb::ConvertExpBetaToMu(double exp_beta, double beta0) {
    return log(1-exp_beta)/(-beta0);
}

double MutationProb::CalculateBeta0(Array4D freq) {

    double beta0 = 1;
    for (size_t i = 0; i < freq.size(); ++i) {
        beta0 -= (freq[i] * freq[i]);
    }
    beta0 = 1.0 / beta0;
    return beta0;

}

