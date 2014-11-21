#include "F81.h"
#include <iostream>
#include "EvolutionModel.h"

EvolutionModel::EvolutionModel(double mu): mu(mu) {}

EvolutionModel::~EvolutionModel() {}

void EvolutionModel::UpdateMu(double mu) {
    this->mu = mu;
    UpdateTransitionMatrix();

}

MutationMatrix EvolutionModel::GetTranstionMatirxAToD() {
    return transition_matrix_a_to_d;
}