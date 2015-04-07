#pragma once
#ifndef _MODEL_H
#define _MODEL_H

//#pragma GCC diagnostic push
//// in reality, you will likely need to disable *more* than Wmultichar
//#pragma GCC diagnostic ignored "-Weverything"

//#pragma GCC diagnostic pop

#include <stdint.h>
#include <vector>

#include "Eigen/Dense"
union ReadData{
    uint64_t key;
    uint16_t reads[4];

};

typedef std::vector<ReadData> ReadDataVector;

struct ModelInput{// Can probably stand to lose this, started out more complex..

    uint16_t reference;
    ReadDataVector all_reads;


    ModelInput() : reference(-1){
    }

    ModelInput(uint read_data_count) : reference(-1){
        all_reads = ReadDataVector(read_data_count, ReadData{0}     );
    }

    ModelInput(uint16_t reference0, ReadDataVector &all_reads0) : reference(reference0), all_reads(all_reads0) {
    }
};


typedef std::vector<ModelInput> GenomeData;

struct ModelParams{
    double theta;               //
    std::vector<double> nuc_freq;    //ACGT
    double mutation_rate;       //
    double error_prob;          // Sequencing error-rate 
    double phi_haploid;         // Overdispersion for haploid sequencing
    double phi_diploid;         // Overdispersion for diploid sequencing
};


typedef Eigen::Array4d HaploidProbs;
typedef Eigen::Array<double, 16, 1> DiploidProbs;
typedef Eigen::Array<double, 16, 4> MutationMatrix;

typedef std::array<double, 10> DiploidProbsIndex10;

HaploidProbs HaploidSequencing(const ModelParams &params, int ref_allele, ReadData data);

DiploidProbs DiploidSequencing(const ModelParams &params, int ref_allele, ReadData data);

MutationMatrix MutationAccumulation(const ModelParams &params, bool and_mut);

DiploidProbs DiploidPopulation(const ModelParams &params, int ref_allele);

void SimulateGenomeData(GenomeData &genome_data, int descendant_count, size_t fake_sample_count, double fake_prop);


double TetMAProbOneMutation(const ModelParams &params, const ModelInput site_data);
double TetMAProbability(const ModelParams &params, const ModelInput site_data);

void SimulateGenomeData(GenomeData &genome_data, int descendant_count, size_t fake_sample_count, double fake_prop = 0);
#endif
