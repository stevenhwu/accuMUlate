//
// Created by steven on 3/31/15.
//
#pragma once
#ifndef _SEQUENCING_FACTORY_H_
#define _SEQUENCING_FACTORY_H_

#include "model.h"
#include "site_genotypes.h"
#include "sequence_prob_v1.h"
#include "site_genotype_index.h"

class SequencingFactory{


public:

    SequencingFactory(ModelParams const &model_params) ;

    void CreateSequenceProbsVector(std::vector<SiteGenotypesIndex> &sp, GenomeData &data);

    const std::vector<HaploidProbs> GetConvertIndexKeyToHaploid();
    const std::vector<DiploidProbsIndex10> GetConvertIndexKeyToDiploidIndex10();

    const std::vector<HaploidProbs> GetConvertIndexKeyToHaploidUnnormalised();
    const std::vector<DiploidProbsIndex10> GetConvertIndexKeyToDiploidIndex10Unnormalised();

//    void CreateSequenceProbsVector(std::vector<SiteGenotypes> &sp, GenomeData &data);

//    void CreateSequenceProbV1(std::vector<SequenceProb> &sp, GenomeData &data);


//    std::vector<DiploidProbs> &GetConvertIndexKeyToDiploid();
//    std::array<DiploidProbs, 4> &GetRefDiploidProbs();


private:

    ModelParams model_params;
    double phi_haploid;
    double phi_diploid;
    double error_prob;
    double theta;
    std::vector<double> frequency_prior;
    Array10D ancestor_prior;

    double haploid_alphas[4][4];


    std::array<DiploidProbs, 4> ref_diploid_probs;

//    std::unordered_map<uint64_t, uint> map_rd_to_index;
//    std::array<std::unordered_map<uint64_t, uint>, 4> map_ancestor_to_index;

    std::vector<HaploidProbs> convert_index_key_to_haploid;
    std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10;

    std::vector<HaploidProbs> convert_index_key_to_haploid_unnormalised;
    std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10_unnormalised;

//
//    void CalculateDescendantGenotypes(SiteGenotypes &seq_prob);
//    void CalculateAncestorGenotype(SiteGenotypes &seq_prob);
//
//    void CalculateDescendantGenotypesIndex(SiteGenotypesIndex &seq_prob);
//    void CalculateAncestorGenotypeIndex(SiteGenotypesIndex &seq_prob);

    DiploidProbs CreateRefDiploidProbs(int ref_allele);
    DiploidProbs DiploidSequencing(ReadData const &data);
    HaploidProbs HaploidSequencing(ReadData const &data);


    void CalculateAncestorPrior();

    DiploidProbsIndex10 ConvertDiploid16ToDiploid10(DiploidProbs probs, int reference);

};
#endif //_ACCUMULATE_SEQUENCING_FACTORY_H_
