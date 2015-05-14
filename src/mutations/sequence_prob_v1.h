/*
 * sequence_prob.h
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#pragma once
#ifndef SEQUENCE_PROB_V1_H_
#define SEQUENCE_PROB_V1_H_

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
#include <mutations/model.h>
#include "mutation_prob.h"
#include "distributions/DirichletMultinomialDistribution.h"
#include "evolution_models/EvolutionModel.h"
#include "lookup.h"
#include <set>
#include <unordered_map>

//
//extern const int ANCESTOR_COUNT;
//extern const int BASE_COUNT;
//
//extern Eigen::IOFormat nice_row;


class SequenceProb {
public:
//    SequenceProb();

    SequenceProb();

    static std::array<DiploidProbs, 4> DiploidPopulationFactory(ModelParams const model_params);
    static HaploidProbs HaploidProbsFactory(ReadData const &data);


    static void CreateSequenceProbV1(std::vector<SequenceProb> &sp, GenomeData &genome_data, ModelParams model_params) {
        sp.reserve(genome_data.size());
        for (size_t i = 0; i < genome_data.size(); ++i) {
            sp.emplace_back(genome_data[i], model_params);
        }
    }

    void SetDescendantGenotypes(std::vector<HaploidProbs> &vector) {
        all_descendant_genotypes = vector;
    }
    void SetAncestorGenotypes(DiploidProbs &anc_genotype){
        ancestor_genotypes = anc_genotype;
    }
    void AddModel(ModelParams const &model_params) {

        phi_haploid = model_params.phi_haploid;
        phi_diploid = model_params.phi_diploid;
        error_prob = model_params.error_prob;
        theta = model_params.theta;
        frequency_prior = model_params.nuc_freq;

    }
//    SequenceProb(){};
    void SetupDiploid(ModelInput const &site_data) {
        size_t i = 0;

        DiploidProbs pop_genotypes = DiploidPopulation(site_data.reference);
        ancestor_data = site_data.all_reads[i];
        ancestor_genotypes = DiploidSequencing(ancestor_data);
        ancestor_genotypes *= pop_genotypes;

        descendant_count = site_data.all_reads.size() - 1;
        descendant_index.assign(descendant_count, -1);

        all_descendant_data.reserve(descendant_count);
        all_descendant_genotypes.reserve(descendant_count);

        for (i = 1; i < (descendant_count + 1); ++i) {
            ReadData data = site_data.all_reads[i];
            all_descendant_data.push_back(data);

        }
    }

    void SetupHaploid() {

    }

//    SequenceProb(const SequenceProb& s);
//    SequenceProb& operator=(const SequenceProb& s);
//    SequenceProb(SequenceProb&& s);

//    SequenceProb(const SequenceProb &s) {
//        std::cout << "Copy Constructor" << std::endl;
//
//    }
//    SequenceProb& operator=(const SequenceProb &s) {
//        std::cout << "Copy assignment operator" << std::endl;
//
//    }
//
//    SequenceProb(SequenceProb &&s) {
//    	std::cout << "Move Constructor" << std::endl;
//    }
//
//    SequenceProb& operator= (SequenceProb&& s) {
//        std::cout << "Move assignment operator" << std::endl;
//    }


    SequenceProb(ModelInput const &site_data, ModelParams const &model_params);
    ~SequenceProb();


//    void UpdateMuProb(MutationProb muProb);
//    void UpdateTransitionMatrix(EvolutionModel evo_model);
    int GetDescendantCount();
    DiploidProbs GetAncestorGenotypes();
    HaploidProbs GetDescendantGenotypes(int descent_index);
    std::vector<HaploidProbs> GetDescendantGenotypes();

    ReadData GetDescendantReadData(int descent_index);
    uint64_t GetDescendantReadDataKey(int descent_index);


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


    ReadData ancestor_data;
    ReadDataVector all_descendant_data;

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

};

#endif /* SEQUENCE_PROB_H_ */
