/*
 * SequenceProb.h
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */

#ifndef SEQUENCEPROB_H_
#define SEQUENCEPROB_H_

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
#include "evolution_models/EvolutionModel.h"
#include "Lookup.h"


extern const int ANCESTOR_COUNT;
extern const int BASE_COUNT;

extern Eigen::IOFormat nice_row;


class SequenceProb {
public:
    static array<DiploidProbs, 4> DiploidPopulationFactory(ModelParams const model_params);


    SequenceProb(ModelInput const site_data, ModelParams const model_params);
    ~SequenceProb();


//    void UpdateMuProb(MutationProb muProb);
//    void UpdateTransitionMatrix(EvolutionModel evo_model);
    int GetDescendantCount();
    DiploidProbs GetAncestorGenotypes();
    HaploidProbs GetDescendantGenotypes(int descent_index);
    std::vector<HaploidProbs> GetDescendantGenotypes();

    void PrintReads(ReadData data);
//    double GetLikelihood();

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

    vector<double> frequency_prior;
    Array10D ancestor_prior;

    int descendant_count;
    double phi_haploid;
    double phi_diploid;
    double error_prob;
    double theta;



};

#endif /* SEQUENCEPROB_H_ */
