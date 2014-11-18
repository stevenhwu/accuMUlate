/*
 * SequencingProb.h
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */

#ifndef SEQUENCINGPROB_H_
#define SEQUENCINGPROB_H_

#include <vector>
#include "model.h"
#include "MutationProb.h"
#include "distributions/DirichletMultinomialDistribution.h"


typedef Eigen::Array<double, 10, 1> Array10D;
typedef HaploidProbs Array4D;



class SequenceProb{
	DiploidProbs anc_genotypes;
	DiploidProbs pop_genotypes;
	std::vector<HaploidProbs> all_hap;
    std::vector<HaploidProbs> all_normalised_hap;
	ModelParams params;

	MutationMatrix non_mutation;
	MutationMatrix mutation;

	double likelihood;
    double mu;
    double beta;
	ModelInput data;

    Array10D frequency_prior;


    int index_convert_16_to_10[4][4] = {
            {0,1,2,3},
            {1,4,5,6},
            {2,5,7,8},
            {3,6,8,9}

    };
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


    Array10D CalculateAncestorToDescendant(Array4D prob_reads_given_descent);
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


};


const double four_third= 4.0/3.0;
const double quarter = 1.0/4.0;
const double three_quarter = 3.0/4.0;
#endif /* SEQUENCINGPROB_H_ */
