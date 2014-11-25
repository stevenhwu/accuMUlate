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

    ReadData ancestor_data;
    ReadDataVector all_descendant_data;

    DiploidProbs ancestor_genotypes;
    std::vector<HaploidProbs> all_descendant_genotypes;

    MutationMatrix transition_matrix_a_to_d;
    MutationRate mutation_rate;

    Array4D frequency_prior;
    Array10D ancestor_prior;

//    std::vector<HaploidProbs> all_normalised_hap;
//    DiploidProbs pop_genotypes;
    double likelihood;
//    MutationMatrix non_mutation;
//    MutationMatrix mutation;
//    double exp_beta;
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


    void CalculateAncestorToDescendant(double &stat_same, double &stat_diff);

    HaploidProbs GetDescendantGenotypes(int descent_index);

    DiploidProbs GetAncestorGenotypes();


    double CalculateExpectedValueForMu(Array10D summary_stat_AtoD);

    double Maximisation(double summery_stat);

    void UpdateTransitionMatrix(EvolutionModel evo_model);

protected:
    DiploidProbs DiploidPopulation(int ref_allele);

    HaploidProbs HaploidSequencing(ReadData data);

    DiploidProbs DiploidSequencing(ReadData data);

//	MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut);
//	MutationMatrix MutationAccumulation2(bool and_mut);


    void PrintReads(ReadData data);

    template<class T>
    T NormaliseLogArray(T result);

private:

    int descendant_count;


    void CalculateDescendantGivenAncestor(int a, HaploidProbs prob_reads_given_descent, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff);
    void CalculateAllDescendantGivenAncestor(int a, double summary_stat_same_ancestor[], double summary_stat_diff_ancestor[], double sum_prob_d[]);


    double phi_haploid;
    double phi_diploid;
    double error_prob;
    double theta;
};

#endif /* SEQUENCINGPROB_H_ */
