#include "JC69.h"

JC69::~JC69() {

}

JC69::JC69(double mu) : EvolutionModel() {

}

void JC69::UpdateMu(double mu) {
    this.mu = mu;

}
