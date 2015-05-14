
#include "JC69.h"

JC69::~JC69() {

}

JC69::JC69(MutationProb mutation_prob) : EvolutionModel(mutation_prob) {
//    cout << "C_JC_PROB\n";
}

JC69::JC69(double mu) : EvolutionModel( MutationProb(mu) ) {
//    cout << "C_JC_mu\n";
////    UpdateTransitionMatrix();

}


void JC69::UpdateTransitionMatrix() {
    UpdateConditionalProbSpecial();
    for (int i : { 0, 1, 2, 3 }) {
        for (int j : { 0, 1, 2, 3 }) {
            int index16 = LookupTable::index_converter_4_4_to_16[i][j];
            int index10 = LookupTable::index_converter_4_4_to_10[i][j];
            for (int k : { 0, 1, 2, 3 }) {
                transition_matrix_a_to_d(index16, k) = conditional_prob_special[index_vector[k][index10]];
            }
        }
    }

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

