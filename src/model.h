#ifndef model_H
#define model_H

#pragma clang diagnostic push
// in reality, you will likely need to disable *more* than Wmultichar
#pragma clang diagnostic ignored "-Wall"
#include "Eigen/Dense"
#pragma clang diagnostic pop


using namespace std;//TODO: maybe remove this later to avoid confusion

union ReadData{
    uint16_t reads[4];
    uint64_t key;
};

typedef vector<ReadData> ReadDataVector;

struct ModelInput{// Can probably stand to lose this, started out more complex..
    uint16_t reference;
    ReadDataVector all_reads;
};

typedef vector<ModelInput> GenomeData;

struct ModelParams{
    double theta;               //
    vector<double> nuc_freq;    //ACGT
    double mutation_rate;       //
    double error_prob;          // Sequencing error-rate 
    double phi_haploid;         // Overdispersion for haploid sequencing
    double phi_diploid;         // Overdispersion for diploid sequencing
};


typedef Eigen::Array4d HaploidProbs;
typedef Eigen::Array<double, 16, 1> DiploidProbs;
typedef Eigen::Array<double, 16, 4> MutationMatrix;



double TetMAProbOneMutation(const ModelParams &params, const ModelInput site_data);
double TetMAProbability(const ModelParams &params, const ModelInput site_data);


#endif
