/*
 * mutation_model.h
 *
 *  Created on: 12/15/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef MUTATION_MODEL_H_
#define MUTATION_MODEL_H_


#include <vector>
#include <iostream>
#include "model.h"
#include "mutation_prob.h"
#include "site_prob.h"
#include "sequence_prob.h"
#include "distributions/DirichletMultinomialDistribution.h"
#include "evolution_models/EvolutionModel.h"
#include "lookup.h"
#include "site_prob.h"

class MutationModel {

public:
    MutationModel(MutationProb mutation_prob, EvolutionModel evo_model);

    virtual ~MutationModel() {
    }

public:
    MutationModel() {
    }

    void CalculateAncestorToDescendant(double &prob_reads, double &all_stats_same, double &all_stats_diff);

    void CalculateAllDescendantGivenAncestor(int a, double &product_prob_given_ancestor, double &summary_stat_same_ancestor, double &summary_stat_diff_ancestor);

    void CalculateOneDescendantGivenAncestor(int anc_index10, HaploidProbs prob_reads_given_descent, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff);

private:

    DiploidProbs ancestor_genotypes;
    std::vector<HaploidProbs> all_descendant_genotypes;

    MutationRate mutation_rate;
    MutationMatrix transition_matrix_a_to_d;

    Array4D frequency_prior;
    Array10D ancestor_prior;

    int descendant_count;


    void UpdateMuProb(MutationProb mutation_prob);

    void UpdateTransitionMatrix(EvolutionModel evo_model);

    void UpdateModel(EvolutionModel &evo_model);


};


#endif //MUTATION_MODEL_H_
