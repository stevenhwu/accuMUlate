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
#include <unordered_map>
#include <map>
#include <mutex>

#include <mutations/model.h>
#include "mutation_prob.h"
#include "site_prob.h"
#include "site_genotypes.h"
#include "distributions/DirichletMultinomialDistribution.h"
//#include "evolution_models/EvolutionModel.h"
#include "lookup.h"
#include "site_prob.h"
#include "sequencing_factory.h"
//#include <boost/functional/hash.hpp>
//#include <boost/unordered_map.hpp>
//#include <sparsehash/dense_hash_map>
//#include <sparsehash/sparse_hash_map>

typedef std::array<std::array<double, 4>, 10> Std2DArray;


class MutationModel {

public:

    MutationModel(EvolutionModel &evo_model0, std::vector<SiteGenotypesIndex> &all_index);

    MutationModel(EvolutionModel &evo_model0);

    MutationModel() {
    }

    virtual ~MutationModel() {
    }

//    void CacheLoopDesAll(int site_index, int anc_index, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);

    void MoveSequenceProb(std::vector<SiteGenotypesIndex> &&all);
//    void MoveSequenceProb(std::vector<SiteGenotypesIndex> all);

    void InitCache();

    void UpdateCache();

    void UpdateOneMinusExpBeta(double oneMinusExpBeta);

    void CalculateAncestorToDescendant(int site_index, double &prob_reads, double &all_stats_diff,
                                       double &log_likelihood_scaler);

    void NoCacheCalculateDes(int site_index, int a, double &product_prob_given_ancestor, double &stat_same,
                             double &summary_stat_diff_ancestor);

    void CalculateOneDescendantGivenAncestor(int anc_index10, HaploidProbs &prob_reads_given_descent,
                                             double &prob_reads_d_given_a, double &summary_stat_same,
                                             double &summary_stat_diff);

    int GetSiteCount() const;
//
////    void UpdateModel(EvolutionModel &evo_model);
////    void CalculateLikelihood(SequenceProbV1 &sequence_prob);
//    void CalculateAllDescendantGivenAncestor(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);

private:

//    std::mutex lock;
//    DiploidProbs ancestor_genotypes;
//    std::vector<HaploidProbs> all_descendant_genotypes;

    EvolutionModel *evo_model;
    double mutation_rate;
    MutationMatrix transition_matrix_a_to_d;

    Array4D frequency_prior;
    Array10D ancestor_prior;
    static std::vector<SiteGenotypesIndex> a;
    static std::vector<SiteGenotypesIndex> all_sequence_prob_index;

    static std::vector<HaploidProbs> convert_index_key_to_haploid;
    static std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10;

    static std::vector<double> convert_index_key_to_haploid_scaler;
    static std::vector<double> convert_index_key_to_diploid_10_scaler;

    static std::vector<HaploidProbs> convert_index_key_to_haploid_unnormalised;
    static std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10_unnormalised;

//    static std::vector<DiploidProbs> convert_index_key_to_diploid;
//    static std::array<DiploidProbs, 4> ref_diploid_probs;



//    std::vector<SiteGenotypes> all_sequence_prob;
//    std::vector<std::array<double, 10>> all_ancestor_genotype_prior;
//    std::unordered_map<uint64_t, HaploidProbs> map_rd_key_to_haploid;



//    std::size_t hash_value(uint16_t read[4]) {
//        std::size_t seed = 0;
//        boost::hash_combine(seed, read[0]);
//        boost::hash_combine(seed, read[1]);
//        boost::hash_combine(seed, read[2]);
//        boost::hash_combine(seed, read[3]);
//        return seed;
//    }

//    std::unordered_map<uint64_t , std::array<std::array<double, 2>, 10> > cache_read_data_to_all; //[key][anc_index]
//    std::array<std::unordered_map<uint64_t, std::array<double, 2>>, 10 > cache_read_data_to_all2; // [anc_index][key]
//    google::sparse_hash_map<uint64_t, std::array<std::array<double, 2>, 10> > cache_read_data_to_all; //[key][anc_index]
//    boost::unordered_map<uint64_t, std::array<std::array<double, 2>, 10> > cache_read_data_to_all_boost; //[key][anc_index]

//    std::unordered_map<uint64_t, int> map_rd_to_index;


//    std::vector<std::array<std::array<double, 2>, 10> > cache_read_data_to_all_index;
//    std::vector<std::array<std::pair<double, double>, 10> > cache_read_data_to_all_index;

    std::array<std::vector<std::pair<double, double>>, 10> cache_read_data_to_all_index_rev;
//    std::array<std::vector<std::array<double, 2>>, 10>  cache_read_data_to_all_index_rev;

//    std::vector<double> summary_stat_diff_vec ;
    int site_count;
    int descendant_count;


public:
    void CacheLoopDesAll2(int anc_index, std::vector<uint32_t> const &aa, double &product_prob_given_ancestor,
                          double &summary_stat_diff_ancestor);

    void CacheLoopDesAll3(int anc_index, std::vector<int> const &aa, double &product_prob_given_ancestor,
                          double &summary_stat_diff_ancestor);

    void CacheLoopDesAll4(int anc_index, std::vector<std::array<std::array<double, 2>, 10> *> &cd,
                          double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);

//    void InitCacheOld1();
//    void SummaryIndexToHaploid();

    static void AddGenotypeFactory(SequencingFactory &factory);

//    void AddSequenceProbOld1(std::vector<SiteGenotypes> &all);
    void CalculateLikelihoodUnnormalised(int site_index, int anc_index, double &product_prob_given_ancestor,
                                         double &summary_stat_diff_ancestor);
};


#endif //MUTATION_MODEL_H_
