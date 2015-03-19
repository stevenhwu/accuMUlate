/*
 * site_prob.h
 *
 *  Created on: Dec 1, 2014
 *      Author: Steven Wu
 */
#pragma once
#ifndef SITE_PROB_H_
#define SITE_PROB_H_

#include <vector>
#include <iostream>
#include "model.h"
#include "mutation_prob.h"
#include "site_prob.h"
#include "sequence_prob.h"
#include "distributions/DirichletMultinomialDistribution.h"
#include "evolution_models/EvolutionModel.h"
#include "lookup.h"

struct ProbTwoStats{

    double prob;
    double stat_same;
    double stat_diff;
};
class SiteProb {


public:
    SiteProb() {
    }

    SiteProb(SequenceProb &sequence_prob, MutationProb &mutation_prob, EvolutionModel &evo_model);

//    SiteProb(ModelInput const site_data, ModelParams const model_params, MutationProb const mutation_prob, EvolutionModel const evo_model);
//    SiteProb(const ModelParams &params, const ModelInput site_data, MutationProb muProb);
//    SiteProb(const ModelInput site_data, const MutationProb muProb, const EvolutionModel evoModel);

    SiteProb(SequenceProb &sequence_prob, EvolutionModel &evo_model);

    ~SiteProb();

    void UpdateMuProb(MutationProb muProb);

    void UpdateTransitionMatrix(EvolutionModel &evo_model);

    void CalculateAncestorToDescendant(double &prob_reads, double &all_stats_same, double &all_stats_diff);

    void CalculateAllDescendantGivenAncestor(int index16, double &sum_prob_d, double &summary_stat_same_ancestor, double &summary_stat_diff_ancestor);

    void CalculateOneDescendantGivenAncestor(int index16, int des_index, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff);


    void UpdateModel(EvolutionModel &evo_model);

protected:

    void PrintReads(ReadData data);


private:

    DiploidProbs ancestor_genotypes;
    std::vector<HaploidProbs> all_descendant_genotypes;
    std::vector<std::array<double, 4>> all_descendant_diff_stats;
    std::array<double, 4> frequency_prior_mutation_rate;

    double mutation_rate;
    MutationMatrix transition_matrix_a_to_d;

    Array4D frequency_prior;
    Array10D ancestor_prior;

    int descendant_count = 0;


};

#endif /* SITE_PROB_H_ */
