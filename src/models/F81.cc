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
        double prob = freqs[i] * (1.0 - exp_beta);
        for (int j : { 0, 1, 2, 3 }) {
            m(i, j) = prob;
        }
        m(i, i) += exp_beta;
    }

    //TODO: use col row to assign them
    conditional_prob = MutationMatrix::Zero();
    for (int i : { 0, 1, 2, 3 }) {
        for (int j : { 0, 1, 2, 3 }) {
            int index16 = i * 4 + j;
            for (int k : { 0, 1, 2, 3 }) {
//                conditional_prob(i * 4 + j, k) = 0.0;
                conditional_prob(index16, k) += 0.5 * m(i, k);
                conditional_prob(index16, k) += 0.5 * m(j, k);
            }
        }
    }
    
    std::cout << conditional_prob << endl;

}


