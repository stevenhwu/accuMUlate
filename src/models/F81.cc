#include <iostream>
#include "F81.h"

F81::~F81() {
}

F81::F81(double mu, vector<double> freq) : EvolutionModel(mu) {
    this->freqs = freq;
    beta0 = 1.0;
	for (auto d : freqs) {
		beta0 -= d * d;
	}
    beta0 = 1/ beta0;
    UpdateConditionalProb();
}


MutationMatrix F81::GetConditionalProb() {
    return conditional_prob;
}

void F81::UpdateConditionalProb() {
    
    exp_beta = exp(-beta0 * mu);

    Eigen::Matrix4d m;
    for (int i : { 0, 1, 2, 3 }) {
        for (int j : { 0, 1, 2, 3 }) {
            m(i, j) = freqs[i] * (1.0 - exp_beta);
        }
        m(i, i) += exp_beta;
    }
//    typedef Eigen::Array<double, 16, 4> MutationMatrix;
    conditional_prob;
    for (int i : { 0, 1, 2, 3 }) {
        for (int j : { 0, 1, 2, 3 }) {
            for (int k : { 0, 1, 2, 3 }) {
                conditional_prob(i * 4 + j, k) = 0.0;
                conditional_prob(i * 4 + j, k) += 0.5 * m(i, k);
                conditional_prob(i * 4 + j, k) += 0.5 * m(j, k);
            }
        }
    }
    std::cout << conditional_prob << endl;

}


