#include <math.h>
#include <iostream>
#include "JC69.h"

JC69::~JC69() {

}

JC69::JC69(double mu) : EvolutionModel(mu) {
    UpdateTransitionMatrix();

}

void JC69::UpdateConditionalProbSpecial(){
    computeExpFourThirdMu();
    conditional_prob_special[0] = probNotEqual();
    conditional_prob_special[1] = probOneEqual();
    conditional_prob_special[2] = probThreeEqual();
}


std::array<double, 3> JC69::GetConditionalProbSpecial() {
    return conditional_prob_special;
}


void JC69::computeExpFourThirdMu(){

    exp_four_third_mu = exp(-FOUR_THIRD*mu);
}

double JC69::probThreeEqual() {
    return QUARTER + THREE_QUARTER* exp_four_third_mu;

}

double JC69::probOneEqual() {
    return QUARTER + QUARTER* exp_four_third_mu;

}


double JC69::probNotEqual() {
    return QUARTER - QUARTER* exp_four_third_mu;
}

const int index_vector[4][10] = { //Same as summary stat now, number of mismatches between D (o) and A(m,f)
        {0, 1, 1, 1, 2, 2, 2, 2, 2, 2},// A
        {2, 1, 2, 2, 0, 1, 1, 2, 2, 2},// C
        {2, 2, 1, 2, 2, 1, 2, 0, 1, 2},// G
        {2, 2, 2, 1, 2, 2, 1, 2, 1, 0} // T
};
const int index_converter_16_to_10[4][4] = {
        {0, 1, 2, 3},
        {1, 4, 5, 6},
        {2, 5, 7, 8},
        {3, 6, 8, 9}

};
void JC69::UpdateTransitionMatrix() {
    UpdateConditionalProbSpecial();
    for (int i : { 0, 1, 2, 3 }) {
        for (int j : { 0, 1, 2, 3 }) {
            int index16 = i * 4 + j;
            int index10 = index_converter_16_to_10[i][j];
            for (int k : { 0, 1, 2, 3 }) {
                transition_matrix_a_to_d(index16, k) = conditional_prob_special[index_vector[k][index10]];
            }
        }
    }


}
