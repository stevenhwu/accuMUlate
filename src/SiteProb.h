/*
 * SiteProb.h
 *
 *  Created on: Dec 1, 2014
 *      Author: Steven Wu
 */

#ifndef SiteProb_H_
#define SiteProb_H_

#include <vector>
#include <iostream>
#include "model.h"
#include "MutationProb.h"
#include "SiteProb.h"
#include "SequenceProb.h"
//#include "distributions/DirichletMultinomialDistribution.h"
#include "evolutionModels/EvolutionModel.h"
#include "Lookup.h"


class SiteProb {


public:
    SiteProb(SequenceProb sequence_prob, MutationProb const mutation_prob, EvolutionModel const evo_model);

//    SiteProb(ModelInput const site_data, ModelParams const model_params, MutationProb const mutation_prob, EvolutionModel const evo_model);
//    SiteProb(const ModelParams &params, const ModelInput site_data, MutationProb muProb);
//    SiteProb(const ModelInput site_data, const MutationProb muProb, const EvolutionModel evoModel);

    ~SiteProb();

    void UpdateMuProb(MutationProb muProb);

    void UpdateTransitionMatrix(EvolutionModel evo_model);


    void CalculateAncestorToDescendant(double &sum_prob, double &stat_same, double &stat_diff);

    void CalculateAllDescendantGivenAncestor(int a, double &sum_prob_d, double &summary_stat_same_ancestor, double &summary_stat_diff_ancestor);

    void CalculateOneDescendantGivenAncestor(int anc_index10, HaploidProbs prob_reads_given_descent, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff);


protected:

    void PrintReads(ReadData data);


private:

    DiploidProbs ancestor_genotypes;
    std::vector<HaploidProbs> all_descendant_genotypes;

    MutationRate mutation_rate;
    MutationMatrix transition_matrix_a_to_d;

    Array4D frequency_prior;
    Array10D ancestor_prior;

    int descendant_count;

};

#endif /* SiteProb_H_ */