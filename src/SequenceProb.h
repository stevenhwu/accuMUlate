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
#include "models/F81.h"

//typedef Eigen::Array<double, 10, 1> Array10D;
//typedef HaploidProbs Array4D;

const int ANCESTOR_COUNT = 10;
const int BASE_COUNT = 4;


const int index_vector[4][10] = { //Same as summary stat now, number of mismatches between D (o) and A(m,f)
        {0, 1, 1, 1, 2, 2, 2, 2, 2, 2},// A
        {2, 1, 2, 2, 0, 1, 1, 2, 2, 2},// C
        {2, 2, 1, 2, 2, 1, 2, 0, 1, 2},// G
        {2, 2, 2, 1, 2, 2, 1, 2, 1, 0} // T
};
const double summary_stat_index_lookup[4][10] ={
        {0, 1, 1, 1, 2, 2, 2, 2, 2, 2},// A
        {2, 1, 2, 2, 0, 1, 1, 2, 2, 2},// C
        {2, 2, 1, 2, 2, 1, 2, 0, 1, 2},// G
        {2, 2, 2, 1, 2, 2, 1, 2, 1, 0} // T

//        {0, 0.5, 0.5, 0.5, 1, 1, 1, 1, 1, 1},// A
//        {1, 0.5, 1, 1, 0, 0.5, 0.5, 1, 1, 1},// C
//        {1, 1, 0.5, 1, 1, 0.5, 1, 0, 0.5, 1},// G
//        {1, 1, 1, 0.5, 1, 1, 0.5, 1, 0.5, 0} // T
};
/*
10 cat version
  AA AC AG AT CC CG CT GG GT TT

  descnt = A
      A  C  G  T
 *AA == !! !! !!
  AC =! =! !! !!
  AG =! !! =! !!
  AT =! !! !! =!
 *CC !! == !! !!
  CG !! =! =! !!
  CT !! =! !! =!
 *GG !! !! == !!
  GT !! !! =! =!
 *TT !! !! !! ==


*/
const int index_converter_16_to_10_single[16] = {
        0,1,2,3,
        1,4,5,6,
        2,5,7,8,
        3,6,8,9,
};
const int index_converter_16_to_10[4][4] = {
        {0,1,2,3},
        {1,4,5,6},
        {2,5,7,8},
        {3,6,8,9}

};
const int index_converter_10_to_16[10] = {//??
        0,1,2, 3,
        5,6, 7,
        10,11,
        15
};
const  string genotype_lookup_10[10]= {
        "AA", "AC", "AG", "AT",
        "CC", "CG", "CT",
        "GG", "GT", "TT"
};

class SequenceProb{



    DiploidProbs anc_genotypes;
	DiploidProbs pop_genotypes;
	std::vector<HaploidProbs> all_hap;
    std::vector<HaploidProbs> all_normalised_hap;
	ModelParams params;

	MutationMatrix non_mutation;
	MutationMatrix mutation;

	double likelihood;
    double beta;
    MutationRate mutation_rate;

	ModelInput data;

    Array4D frequency_prior;
    Array10D ancestor_prior;



//    const char ACGT[] {'A', 'C', 'G', 'T'};
public:
	SequenceProb(){};
	SequenceProb(const ModelParams &params, const ModelInput site_data, MutationProb muProb);
	~SequenceProb();

//	void UpdateMu();
//	void UpdateMu(double mu);
	void UpdateLikelihood();
	void UpdateMuProb(MutationProb muProb);

	double GetLikelihood();

	ModelInput GetData();


    void CountReadToGenotype();


    void CalculateAncestorToDescendant(MutationMatrix &conditional_prob);
    HaploidProbs GetDescendantToReads(int descent_index);

    DiploidProbs GetAncGenotypesToReads();

    static double probNotEqual(double emu);
    static double probOneEqual(double emu);
    static double probThreeEqual(double emu);
    static double computeExpFourThirdMu(double mu);

    double CalculateExpectedValueForMu(Array10D summary_stat_AtoD);

    double Maximisation(double summery_stat);

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

    size_t descendant_count;

    void CalculateDescendantGivenAncestor(MutationMatrix &prob_matrix_a_d, double mu, int a, HaploidProbs prob_reads_given_descent, double &sum_over_probs, double &summary_stat);

    void CalculateAllDescendantGivenAncestor(MutationMatrix &conditional_prob, int a, double summary_stat_d[], double summary_stat_same_ancestor[], double sum_prob_d[]);


};

#endif /* SEQUENCINGPROB_H_ */
