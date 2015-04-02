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
    static void CreateSequenceProb(SequenceProb &sp, ModelInput const &data, ModelParams const params);

    SequencingFactory(ModelParams const &model_params) ;

    SequenceProb CreateSequenceProb(ModelInput const &site_data);


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


    std::array<DiploidProbs, 4> ref_diploid_probs;
    std::vector<HaploidProbs> convert_index_key_to_haploid;
    std::vector<DiploidProbs> convert_index_key_to_diploid;


    void CalculateDescendantGenotypes(SequenceProb &seq_prob);
    void CalculateAncestorGenotype(SequenceProb &seq_prob);


    DiploidProbs CreateRefDiploidProbs(int ref_allele);
    DiploidProbs DiploidSequencing(ReadData const &data);
    HaploidProbs HaploidSequencing(ReadData const &data);


};
#endif //_ACCUMULATE_SEQUENCING_FACTORY_H_
