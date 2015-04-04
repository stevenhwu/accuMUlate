/*
 * sequence_prob.h
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#pragma once
#ifndef SEQUENCE_PROB_H_
#define SEQUENCE_PROB_H_

//#ifdef __GNUC__
//#define DEPRECATED(func) func __attribute__ ((deprecated))
//#elif defined(_MSC_VER)
//#define DEPRECATED(func) __declspec(deprecated) func
//#else
//#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
//#define DEPRECATED(func) func
//#endif


#include <vector>
#include <iostream>
#include <set>
#include <unordered_map>

#include "model.h"
#include "mutation_prob.h"
#include "distributions/DirichletMultinomialDistribution.h"
#include "evolution_models/EvolutionModel.h"
#include "lookup.h"
#include "constant.h"



class SiteGenotypes {
public:

    static
    void SiteGenotypes_ASSERT_GENOTYPES(HaploidProbs expected, HaploidProbs data){
        for (int i = 0; i < expected.size(); ++i) {
            assert( abs(expected[i]- data[i]) < ERROR_THRESHOLD );
            assert( expected[i] == data[i] );
        }
    }
    static void SiteGenotypes_ASSERT_GENOTYPES(DiploidProbs expected, DiploidProbs data){
        for (int i = 0; i < expected.size(); ++i) {
            assert( abs(expected[i]- data[i]) < ERROR_THRESHOLD );
            assert( expected[i] == data[i] );
        }
    }

//    SiteGenotypes();

    static std::array<DiploidProbs, 4> DiploidPopulationFactory(ModelParams const model_params);
    static HaploidProbs HaploidProbsFactory(ReadData const &data);
    static void printReadData(ReadData read_data);



    void AddModel(ModelParams const &model_params) {

//        phi_haploid = model_params.phi_haploid;
        phi_diploid = model_params.phi_diploid;
        error_prob = model_params.error_prob;
        theta = model_params.theta;
        frequency_prior = model_params.nuc_freq;

    }
    SiteGenotypes(){}
    uint16_t reference;
    SiteGenotypes(ModelInput const &site_data) {

        size_t i = 0;
        reference = site_data.reference;

        auto rd_vector = site_data.all_reads;
        ancestor_readdata = rd_vector[i];

        descendant_count = site_data.all_reads.size() - 1;
        descendant_index.assign(descendant_count, -1);
        all_descendant_readdata.reserve(descendant_count);
        all_descendant_genotypes.reserve(descendant_count);

        for (i = 1; i < (descendant_count + 1); ++i) {
            all_descendant_readdata.emplace_back(rd_vector[i]);
        }
    }

    void SetupDiploid(ModelInput const &site_data) {
        size_t i = 0;

        auto rd_vector = site_data.all_reads;
        ancestor_readdata = rd_vector[i];

        DiploidProbs pop_genotypes = DiploidPopulation(site_data.reference);
        ancestor_genotypes = DiploidSequencing(ancestor_readdata);
        ancestor_genotypes *= pop_genotypes;

        descendant_count = site_data.all_reads.size() - 1;
        descendant_index.assign(descendant_count, -1);

        all_descendant_readdata.reserve(descendant_count);
        all_descendant_genotypes.reserve(descendant_count);

        for (i = 1; i < (descendant_count + 1); ++i) {
            all_descendant_readdata.emplace_back(rd_vector[i]);
        }

    }
    uint16_t GetReference(){
        return reference;
    }
    void SetDescendantGenotypes(std::vector<HaploidProbs> &vector) {
        all_descendant_genotypes = vector;
    }
    void SetAncestorGenotypes(DiploidProbs &anc_genotype){
        ancestor_genotypes = anc_genotype;
    }
    ReadData GetAncestorReadData(){
        return ancestor_readdata;
    }
    void SetAncestorIndex(int index){
        ancestor_index = index;
    }
    int GetAncestorIndex(){
        return ancestor_index;
    }
//    SiteGenotypes(const SiteGenotypes& s);
//    SiteGenotypes& operator=(const SiteGenotypes& s);
//    SiteGenotypes(SiteGenotypes&& s);

//    SiteGenotypes(const SiteGenotypes &s) {
//        std::cout << "Copy Constructor" << std::endl;
//
//    }
//    SiteGenotypes & operator=(const SiteGenotypes &s) {
//        std::cout << "Copy assignment operator" << std::endl;
//
//    }
//
//    SiteGenotypes(SiteGenotypes &&s) {
//    	std::cout << "Move Constructor" << std::endl;
//    }
//
//    SiteGenotypes & operator= (SiteGenotypes && s) {
//        std::cout << "Move assignment operator" << std::endl;
//    }


    SiteGenotypes(ModelInput const &site_data, ModelParams const &model_params);
    ~SiteGenotypes();


//    void UpdateMuProb(MutationProb muProb);
//    void UpdateTransitionMatrix(EvolutionModel evo_model);
    int GetDescendantCount();
    DiploidProbs GetAncestorGenotypes();
    HaploidProbs GetDescendantGenotypes(int descent_index);
    std::vector<HaploidProbs> GetDescendantGenotypes();

    ReadData GetDescendantReadData(int descent_index);
    uint64_t GetDescendantReadDataKey(int descent_index);

    void PrintReads(ReadData data);
//    double GetLikelihood();

    ReadDataVector GetDescendantReadData();
    ReadDataVector const & GetDescendantReadData2();
    ReadDataVector const * GetDescendantReadData3();
    void GetDescendantReadDataCOPY(ReadDataVector &all_descendant_data);

    std::unordered_map<uint64_t, int> temp_map;
    std::vector<std::pair<uint64_t, int>> condense_genotype;

    void SetDescendantIndex(int des, int index);
    int GetDescendantIndex(int des);


    const std::vector<int>& GetDescendantIndex();

    void SortIndex();

//protected:
    DiploidProbs DiploidPopulation(int ref_allele);

    HaploidProbs HaploidSequencing(ReadData const &data);

    DiploidProbs DiploidSequencing(ReadData const &data);


    template<class T>
    T ScaleLogArray(T result);


private:


    ReadData ancestor_readdata;
    ReadDataVector all_descendant_readdata;

    DiploidProbs ancestor_genotypes;
    std::vector<HaploidProbs> all_descendant_genotypes;

    std::vector<double> frequency_prior;
    Array10D ancestor_prior;

    size_t descendant_count;
    double phi_haploid;
    double phi_diploid;
    double error_prob;
    double theta;

    std::vector<int> descendant_index;
    int ancestor_index;
};

#endif /* SEQUENCE_PROB_H_ */
