#include "EvolutionModel.h"

EvolutionModel::EvolutionModel(double mu): mu(mu) {

}

EvolutionModel::~EvolutionModel() {}

void EvolutionModel::UpdateMu(double mu) {
    this->mu = mu;
    UpdateConditionalProb();

}


