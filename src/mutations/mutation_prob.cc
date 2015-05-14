/*
 * mutation_prob.cpp
 *
 *  Created on: Nov 8, 2014
 *      Author: Steven Wu
 */

#include "mutation_prob.h"

MutationProb::MutationProb(const ModelParams &model_params) {

    for (int i = 0; i < 4; ++i) {
        frequency_prior[i] = model_params.nuc_freq[i];
    }

    CalculateAncestorPrior();
    beta0 = CalculateBeta0(frequency_prior);
    UpdateMu(model_params.mutation_rate);
}

MutationProb::MutationProb() :MutationProb(1){};
MutationProb::MutationProb(double mu0) {

    for (int i = 0; i < 4; ++i) {
        frequency_prior[i] = 0.25;
    }
    CalculateAncestorPrior();
    beta0 = CalculateBeta0(frequency_prior);
    UpdateMu(mu0);

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
    this->mu = mu0;
    CalculateExpBeta();

}



double MutationProb::CalculateExpBeta(){
    exp_beta = CalculateExpBeta(mu, beta0);
    mutation_rate = 1 - exp_beta;
    return exp_beta;
}

double MutationProb::ConvertOneMinusExpBetaToMu(double one_minus_exp_beta) {
    return ConvertOneMinusExpBetaToMu(one_minus_exp_beta, beta0);
}

double MutationProb::GetMutationRate() {
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

double MutationProb::GetMu() {
    return mu;
}

double MutationProb::GetBeta0() {
    return beta0;
}


//////////////////////////////////////////////////
//// Static methods
//////////////////////////////////////////////////
double MutationProb::CalculateExpBeta(double mu0, double beta0){
    double exp_beta = exp(-beta0 * mu0);
    return exp_beta;
}


double MutationProb::ConvertOneMinusExpBetaToMu(double one_minus_exp_beta, double beta0) {
    return MutationProb::ConvertExpBetaToMu(1-one_minus_exp_beta, beta0);
}

double MutationProb::ConvertExpBetaToMu(double exp_beta0, double beta0) {
    return log(exp_beta0)/(-beta0);
}


double MutationProb::CalculateBeta0(Array4D &freq) {

    double beta0 = 1;
    for (size_t i = 0; i < freq.size(); ++i) {
        beta0 -= (freq[i] * freq[i]);
    }
    beta0 = 1.0 / beta0;
    return beta0;

}
