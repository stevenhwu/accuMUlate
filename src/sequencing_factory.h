//
// Created by steven on 3/31/15.
//

#ifndef _SEQUENCING_FACTORY_H_
#define _SEQUENCING_FACTORY_H_

#include "model.h"
#include "sequence_prob.h"
#include "sequence_prob_v1.h"

class SequencingFactory{


public:


    SequencingFactory(ModelParams const &model_params) ;

    SequenceProb CreateSequenceProb(ModelInput const &site_data);

    SequenceProb & InitSequenceProb(SequenceProb &seq_prob, ModelInput &site_data);

    void CreateSequenceProbsVector(std::vector<SequenceProb> &sp, GenomeData &data);
    void CreateSequenceProbV1(std::vector<SequenceProbV1> &sp, GenomeData &data);

private:

    ModelParams model_params;
    double phi_haploid;
    double phi_diploid;
    double error_prob;
    double theta;
    std::vector<double> frequency_prior;

    double haploid_alphas[4][4];


//    ReadData &ancestor_data;
//    DiploidProbs ancestor_genotypes;
    unsigned long descendant_count;

//    ReadDataVector all_descendant_data;
//    std::vector<HaploidProbs> all_descendant_genotypes;
//    Array10D ancestor_prior;
//    std::vector<int> descendant_index;
    std::vector<HaploidProbs> convert_index_key_to_haploid;
    std::unordered_map<uint64_t, HaploidProbs> map_rd_key_to_haploid;//remove later
    void CalculateHaploidProb(SequenceProb &seq_prob, ModelInput &site_data);

    HaploidProbs HaploidSequencing(ReadData const &data);
};
#endif //_ACCUMULATE_SEQUENCING_FACTORY_H_
