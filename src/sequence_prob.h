/*
 * sequence_prob.h
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#pragma once
#ifndef SEQUENCE_PROB_H_
#define SEQUENCE_PROB_H_

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED(func) func
#endif


#include <vector>
#include <iostream>
#include "model.h"
#include "mutation_prob.h"
#include "distributions/DirichletMultinomialDistribution.h"
#include "evolution_models/EvolutionModel.h"
#include "lookup.h"
#include <set>
#include <unordered_map>

extern const int ANCESTOR_COUNT;
extern const int BASE_COUNT;

extern Eigen::IOFormat nice_row;


class SequenceProb {
public:
    static std::array<DiploidProbs, 4> DiploidPopulationFactory(ModelParams const model_params);
    static void printReadData(ReadData read_data);

    SequenceProb(ModelInput const site_data, ModelParams const model_params);
    ~SequenceProb();


//    void UpdateMuProb(MutationProb muProb);
//    void UpdateTransitionMatrix(EvolutionModel evo_model);
    int GetDescendantCount();
    DiploidProbs GetAncestorGenotypes();
    HaploidProbs GetDescendantGenotypes(int descent_index);
    std::vector<HaploidProbs> GetDescendantGenotypes();

    ReadData GetDescendantReadData(int descent_index);

    void PrintReads(ReadData data);
//    double GetLikelihood();

    ReadDataVector GetDescendantReadData();
    ReadDataVector const & GetDescendantReadData2();
    ReadDataVector const * GetDescendantReadData3();
    void GetDescendantReadDataCOPY(ReadDataVector &all_descendant_data);

    std::unordered_map<uint64_t, int> temp_map;
    std::vector<std::pair<uint64_t, int>> condense_genotype;
protected:
    DiploidProbs DiploidPopulation(int ref_allele);

    HaploidProbs HaploidSequencing(ReadData data);

    DiploidProbs DiploidSequencing(ReadData data);



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


};

#endif /* SEQUENCE_PROB_H_ */
