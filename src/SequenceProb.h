/*
 * SequencingProb.h
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */

#ifndef SEQUENCINGPROB_H_
#define SEQUENCINGPROB_H_

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
#include "MutationProb.h"
#include "distributions/DirichletMultinomialDistribution.h"
#include "models/EvolutionModel.h"
#include "Lookup.h"
//typedef Eigen::Array<double, 10, 1> Array10D;
//typedef HaploidProbs Array4D;

class SequenceProb {


    DiploidProbs anc_genotypes;
    DiploidProbs pop_genotypes;
    std::vector<HaploidProbs> all_hap;
    std::vector<HaploidProbs> all_normalised_hap;
    ModelParams params;

    MutationMatrix non_mutation;
    MutationMatrix mutation;

    MutationMatrix transition_matrix_a_to_d;
    double likelihood;
    double beta;
    MutationRate mutation_rate;

    ModelInput data;

    Array4D frequency_prior;
    Array10D ancestor_prior;



//    c char ACGT[] {'A', 'C', 'G', 'T'};
public:
    SequenceProb() {
    };

    SequenceProb(const ModelParams &params, const ModelInput site_data, MutationProb muProb);

    ~SequenceProb();

//	void UpdateMu();
//	void UpdateMu(double mu);
    void UpdateLikelihood();

    void UpdateMuProb(MutationProb muProb);

    double GetLikelihood();

    ModelInput GetData();


    void CountReadToGenotype();


    void CalculateAncestorToDescendant();

    HaploidProbs GetDescendantToReads(int descent_index);

    DiploidProbs GetAncGenotypesToReads();


    double CalculateExpectedValueForMu(Array10D summary_stat_AtoD);

    double Maximisation(double summery_stat);

    void UpdateTransitionMatrix(EvolutionModel f81);

protected:
    DiploidProbs DiploidPopulation(int ref_allele);

    HaploidProbs HaploidSequencing(ReadData data);

    DiploidProbs DiploidSequencing(ReadData data);

//	MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut);
//	MutationMatrix MutationAccumulation2(bool and_mut);




    ReadData ancestor;
//    vector<ReadData> all_descendant;
    ReadDataVector all_descendant;

    void PrintReads(ReadData data);

    template<class T>
    T NormaliseLogArray(T result);

private:

    int descendant_count;


    void CalculateDescendantGivenAncestor(int a, HaploidProbs prob_reads_given_descent, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff);
    void CalculateAllDescendantGivenAncestor(int a, double summary_stat_same_ancestor[], double summary_stat_diff_ancestor[], double sum_prob_d[]);


};

#endif /* SEQUENCINGPROB_H_ */
