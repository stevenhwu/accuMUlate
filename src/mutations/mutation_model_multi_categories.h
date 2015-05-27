/*
 * mutation_model_multi_categories.cc
 *
 *  Created on: 4/10/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/10/15.
//

#pragma once
#ifndef _ACCUMULATE_MUTATION_MODEL_MULTI_CATEGORIES_H_
#define _ACCUMULATE_MUTATION_MODEL_MULTI_CATEGORIES_H_


#include <vector>
#include <iostream>
#include <unordered_map>
#include <map>
//#include <evolution_models/F81.h>

#include <mutations/model.h>
#include "mutation_prob.h"
#include "site_prob.h"
#include "site_genotypes.h"
#include "distributions/DirichletMultinomialDistribution.h"
#include "evolution_models/EvolutionModel.h"
#include "lookup.h"
#include "site_prob.h"
#include "sequencing_factory.h"

class MutationModelMultiCategories {


public:

    MutationModelMultiCategories(int num_categories, SequencingFactory &sf);

    MutationModelMultiCategories(int num_categories, EvolutionModel  &evo_model0);

    MutationModelMultiCategories() {
    }

    MutationModelMultiCategories(SequencingFactory &factory);


    MutationModelMultiCategories(int num_categories, EvolutionModel &evo_model0, SequencingFactory &factory);

    virtual ~MutationModelMultiCategories() {
    }
//TODO: Current plan, one EvolutionModel,  one sets of data. Two/multi updateExpBeta => Multi mutation rate, multi transition matrix
    void AddEvolutionModels(EvolutionModel  &evo_model0);

    void MoveSequenceProb(std::vector<SiteGenotypesIndex> &&all);

    void AddGenotypeFactory(SequencingFactory &factory);

    void InitCache(int category_index);

    void UpdateCache(int category_index);

    void UpdateOneMinusExpBeta(int category_index, double oneMinusExpBeta);

    void CalculateAncestorToDescendant(int category_index, int site_index, double &prob_reads,
                                                                     double &all_stats_diff, double &log_likelihood_scaler);

    void NoCacheCalculateDes(int categories_index, int site_index, int anc_index10,
                                                           double &product_prob_given_ancestor, double &stat_same,
                                                           double &summary_stat_diff_ancestor);

    void CalculateOneDescendantGivenAncestor(int category_index, int anc_index10,
                                                                           HaploidProbs &prob_reads_given_descent,
                                                                           double &prob_reads_d_given_a,
                                                                           double &summary_stat_same,
                                                                           double &summary_stat_diff);

    int GetCategoriesCount() const;
    int GetSiteCount() const;
//
////    void UpdateModel(EvolutionModel &evo_model);
////    void CalculateLikelihood(SequenceProbV1 &sequence_prob);
//    void CalculateAllDescendantGivenAncestor(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);

    int GetDescendantCount(int site_index);

private:

//    DiploidProbs ancestor_genotypes;
//    std::vector<HaploidProbs> all_descendant_genotypes;

    EvolutionModel *evo_model;

    std::vector<double> mutation_rate_single;
    std::vector<MutationMatrix> transition_matrix_a_to_d_single;

//    std::vector<double> mutation_rate;
//    std::vector<MutationMatrix> transition_matrix_a_to_d;

    Array4D frequency_prior;
    Array10D ancestor_prior;

    std::vector<SiteGenotypesIndex> all_sequence_prob_index;

    std::vector<HaploidProbs> convert_index_key_to_haploid;
    std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10;

    std::vector<double> convert_index_key_to_haploid_scaler;
    std::vector<double> convert_index_key_to_diploid_10_scaler;

    std::vector<HaploidProbs> convert_index_key_to_haploid_unnormalised;
    std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10_unnormalised;
    std::vector<std::array<std::vector<std::pair<double, double>>, 10>> cache_read_data_to_all_index_rev;//TODO: turn this into a class

    int site_count;
//    int descendant_count;



    void CacheLoopDesAll2(int category_index, int anc_index, const std::vector<uint32_t> &aa,
                                                        double &product_prob_given_ancestor,
                                                        double &summary_stat_diff_ancestor);


    void CalculateLikelihoodUnnormalised(int site_index, int anc_index, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);

    int categories_count;
};


#endif //_ACCUMULATE_MUTATION_MODEL_MULTI_CATEGORIES_H_
