#include <iostream>
#include "EvolutionModel.h"

//EvolutionModel::EvolutionModel(double mu):  mu_prob(MutationProb(mu)) {
////    mu_prob = MutationProb(mu);
//    cout << "C_mu:\t" << mu << "\t" << mu_prob.GetMu() << endl;
//
//}

EvolutionModel::EvolutionModel(MutationProb mu_prob) : mu_prob(mu_prob) {

    this->mu = mu_prob.GetMu();
}

EvolutionModel::~EvolutionModel() {}



void EvolutionModel::UpdateMu(double mu0) {
    mu_prob.UpdateMu(mu0);
    this->mu = mu0;
//    cout << "new mu:" << mu << endl;
    UpdateTransitionMatrix();

}

MutationMatrix EvolutionModel::GetTranstionMatirxAToD() {
    return transition_matrix_a_to_d;
}

double EvolutionModel::GetMutationRate() {
//    cout << mu_prob.beta0 << endl;
//    cout << mu_prob.GetExpBeta() << endl;
    return mu_prob.GetMutationRate();

}

MutationProb EvolutionModel::GetMutationProb() {
    return mu_prob;
}

MutationProb EvolutionModel::GetMutationProb() const {
    MutationProb m2 = mu_prob;
    return m2;
}

