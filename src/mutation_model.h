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

#include "model.h"
#include "mutation_prob.h"
#include "site_prob.h"
#include "sequence_prob.h"
#include "distributions/DirichletMultinomialDistribution.h"
#include "evolution_models/EvolutionModel.h"
#include "lookup.h"
#include "site_prob.h"
//#include <boost/functional/hash.hpp>
//#include <boost/unordered_map.hpp>
//#include <sparsehash/dense_hash_map>
//#include <sparsehash/sparse_hash_map>

typedef std::array<std::array<double, 4>, 10> Std2DArray;



class MutationModel {

public:

    MutationModel(EvolutionModel &evo_model0);
    MutationModel() {
    }
    virtual ~MutationModel() {
    }

    void CacheLoopDesAll(int site_index, int anc_index, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);

    void AddSequenceProb(std::vector<SequenceProb> &all);

    void InitCache();
    void UpdateCache();

    void UpdateExpBeta(double expBeta);


    void CalculateAncestorToDescendant(int site_index, double &prob_reads, double &all_stats_diff);

    void NoCacheCalculateDes(int site_index, int a, double &product_prob_given_ancestor, double &stat_same, double &summary_stat_diff_ancestor);
    void CalculateOneDescendantGivenAncestor(int anc_index10, HaploidProbs &prob_reads_given_descent, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff);

    int GetSiteCount() const;

//    void UpdateModel(EvolutionModel &evo_model);
//    void CalculateLikelihood(SequenceProb &sequence_prob);
//    void CalculateAllDescendantGivenAncestor(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);

private:

//    DiploidProbs ancestor_genotypes;
//    std::vector<HaploidProbs> all_descendant_genotypes;

    EvolutionModel *evo_model;
    double mutation_rate;
    MutationMatrix transition_matrix_a_to_d;

    Array4D frequency_prior;
    Array10D ancestor_prior;

    std::vector<SequenceProb> all_sequence_prob;
    std::vector<std::array<double, 10>> all_ancestor_genotype_prior;
    std::unordered_map<uint64_t, HaploidProbs> map_rd_key_to_haploid;




    std::unordered_map<uint64_t , std::array<std::array<double, 2>, 10> > cache_read_data_to_all; //[key][anc_index]
//    std::array<std::unordered_map<uint64_t, std::array<double, 2>>, 10 > cache_read_data_to_all2; // [anc_index][key]
//    google::sparse_hash_map<uint64_t, std::array<std::array<double, 2>, 10> > cache_read_data_to_all; //[key][anc_index]
//    boost::unordered_map<uint64_t, std::array<std::array<double, 2>, 10> > cache_read_data_to_all_boost; //[key][anc_index]

    std::unordered_map<uint64_t, int> map_rd_to_index;
    std::vector<HaploidProbs> map_index_key_to_haploid;


//    std::vector<std::array<std::array<double, 2>, 10> > cache_read_data_to_all_index;
    std::vector<std::array<std::pair<double, double>, 10> > cache_read_data_to_all_index;

    std::array<std::vector<std::pair<double, double>>, 10>  cache_read_data_to_all_index_rev;
//    std::array<std::vector<std::array<double, 2>>, 10>  cache_read_data_to_all_index_rev;

    std::vector<double> summary_stat_diff_vec ;
    int site_count;
    int descendant_count;



public:

//    void UpdateMuProb(MutationProb mutation_prob);
//
//    void UpdateTransitionMatrix(EvolutionModel evo_model);

//    std::unordered_map<double, std::array<double, 4>> cache_data;
//    std::unordered_map<double,  Std2DArray> cache_data_transition;

//    ReadDataVector all_descendant_data;

//    std::unordered_map<uint64_t, int> map;
//    std::vector<std::pair<uint64_t, int>> vec;


//    void CalculateAncestorToDescendantOrig(double &prob_reads, double &all_stats_same, double &all_stats_diff);
//
//    void CalculateAllDescendantGivenAncestorOrig(int a, double &product_prob_given_ancestor, double &summary_stat_same_ancestor, double &summary_stat_diff_ancestor);
//
//    void CalculateOneDescendantGivenAncestorOrig(int anc_index10, HaploidProbs prob_reads_given_descent, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff);
//
//    void CalculateLikelihoodOriginal(SequenceProb &sequence_prob, double &prob, double &stat_same, double &stat_diff);
//
//    void CalculateOneDescendantGivenAncestorCache(int anc_index10, HaploidProbs prob_reads_given_descent, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff);

//    void CacheLoopDes(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);

//    void CacheLoopMap(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);

//    void AddCache2(std::unordered_map<uint64_t, std::array<std::array<double, 2>, 10> > &all_2);
//            std::unordered_map<uint64_t, std::array<double, 10>> prob, std::unordered_map<uint64_t, std::array<double, 10>> stat,
//            std::unordered_map<uint64_t, std::pair<std::array<double, 10>, std::array<double, 10> > > all,

//    void AddCache(std::unordered_map<double, std::array<double, 4>> &allDMaps4,
//            std::unordered_map<double,  Std2DArray>  &cache_data_transition);

    void CacheLoopDesAll2(int anc_index, std::vector<int> const &aa, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);
    void CacheLoopDesAll3(int anc_index, std::vector<int> const &aa, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor);
    void CacheLoopDesAll4(int anc_index, std::vector<std::array<std::array<double, 2>, 10>*> &cd, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) ;
};


#endif //MUTATION_MODEL_H_
