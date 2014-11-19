#include <math.h>
#include <iostream>
#include "JC69.h"

JC69::~JC69() {

}

JC69::JC69(double mu) : EvolutionModel(mu) {
    UpdateConditionalProb();
}

void JC69::UpdateConditionalProb(){
    computeExpFourThirdMu();
    conditional_prob[0] = probNotEqual();
    conditional_prob[1] = probOneEqual();
    conditional_prob[2] = probThreeEqual();
}



//double JC69::GetConditionalProb(int i) {
//    return conditional_prob[i];
//}

std::array<double, 3> JC69::GetConditionalProbSpecial() {
    return conditional_prob;
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
