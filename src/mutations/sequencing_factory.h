//
// Created by steven on 3/31/15.
//
#pragma once
#ifndef _SEQUENCING_FACTORY_H_
#define _SEQUENCING_FACTORY_H_

#include <mutations/model.h>
#include "site_genotypes.h"
#include "sequence_prob_v1.h"
#include "site_genotype_index.h"

class SequencingFactory{


public:

    SequencingFactory(ModelParams const &model_params) ;

    void CreateSequenceProbsVector(std::vector<SiteGenotypesIndex> &sgi, GenomeData &data);
    void CreateSequenceProbsVector(std::vector<SiteGenotypesIndex> &sgi, ModelInput &data);
    void CreateSequenceProbsVector(GenomeData &data);

    const std::vector<HaploidProbs> RemoveConvertIndexKeyToHaploid();
    const std::vector<DiploidProbsIndex10> RemoveConvertIndexKeyToDiploidIndex10();

    std::vector<double> && RemoveConvertIndexKeyToDiploidIndex10Scaler();
    std::vector<double> && RemoveConvertIndexKeyToHaploidScaler();

    std::vector<SiteGenotypesIndex> &&RemoveSiteGenotypeIndexVector();

//    const std::vector<HaploidProbs> RemoveConvertIndexKeyToHaploidUnnormalised();
//    const std::vector<DiploidProbsIndex10> RemoveConvertIndexKeyToDiploidIndex10Unnormalised();

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

    std::vector<SiteGenotypesIndex> sgi;
    std::array<DiploidProbs, 4> ref_diploid_probs;

    std::vector<HaploidProbs> convert_index_key_to_haploid;
    std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10;

    std::vector<double> convert_index_key_to_haploid_scaler;
    std::vector<double> convert_index_key_to_diploid_10_scaler;

//    std::vector<HaploidProbs> convert_index_key_to_haploid_unnormalised;
//    std::vector<DiploidProbsIndex10> convert_index_key_to_diploid_10_unnormalised;


    std::unordered_map<uint64_t, uint32_t> map_rd_to_index;
    std::array<std::unordered_map<uint64_t, uint32_t>, 4> map_ancestor_to_index;

    std::unordered_map<int, int> map_des_count;
    //
//    void CalculateDescendantGenotypes(SiteGenotypes &seq_prob);
//    void CalculateAncestorGenotype(SiteGenotypes &seq_prob);
//
//    void CalculateDescendantGenotypesIndex(SiteGenotypesIndex &seq_prob);
//    void CalculateAncestorGenotypeIndex(SiteGenotypesIndex &seq_prob);
    void CalculateAncestorPrior();

    DiploidProbs CreateRefDiploidProbs(int ref_allele);
    DiploidProbs DiploidSequencing(ReadData const &data);
    HaploidProbs HaploidSequencing(ReadData const &data);

    DiploidProbsIndex10 ConvertDiploid16ToDiploid10(DiploidProbs probs, int reference);

};
#endif //_ACCUMULATE_SEQUENCING_FACTORY_H_
