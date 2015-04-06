//
// Created by steven on 3/31/15.
//

#ifndef _SEQUENCING_FACTORY_H_
#define _SEQUENCING_FACTORY_H_

#include "model.h"
#include "site_genotypes.h"
#include "sequence_prob_v1.h"

class SequencingFactory{


public:
    static void CreateSequenceProb(SiteGenotypes &sp, ModelInput const &data, ModelParams const params);

    SequencingFactory(ModelParams const &model_params) ;

    SiteGenotypes CreateSequenceProb(ModelInput const &site_data);


    void CreateSequenceProbsVector(std::vector<SiteGenotypes> &sp, GenomeData &data);
    void CreateSequenceProbV1(std::vector<SequenceProb> &sp, GenomeData &data);

    std::vector<HaploidProbs> &GetConvertIndexKeyToHaploid();
    std::vector<DiploidProbs> &GetConvertIndexKeyToDiploid();
    std::array<DiploidProbs, 4> &GetRefDiploidProbs();
private:

    ModelParams model_params;
    double phi_haploid;
    double phi_diploid;
    double error_prob;
    double theta;
    std::vector<double> frequency_prior;

    double haploid_alphas[4][4];


    std::array<DiploidProbs, 4> ref_diploid_probs;


    std::unordered_map<uint64_t, int> map_rd_to_index;
    std::unordered_map<uint64_t, int> map_ancestor_to_index;

    std::vector<HaploidProbs> convert_index_key_to_haploid;
    std::vector<DiploidProbs> convert_index_key_to_diploid;


    void CalculateDescendantGenotypes(SiteGenotypes &seq_prob);
    void CalculateAncestorGenotype(SiteGenotypes &seq_prob);


    DiploidProbs CreateRefDiploidProbs(int ref_allele);
    DiploidProbs DiploidSequencing(ReadData const &data);
    HaploidProbs HaploidSequencing(ReadData const &data);


};
#endif //_ACCUMULATE_SEQUENCING_FACTORY_H_
