#pragma once
#ifndef _MODEL_H
#define _MODEL_H

//#pragma GCC diagnostic push
//// in reality, you will likely need to disable *more* than Wmultichar
//#pragma GCC diagnostic ignored "-Weverything"

//#pragma GCC diagnostic pop

//#include <stdint.h>
#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <map>
#include <fstream>
#include <memory>
#include <random>
#include <Eigen/Dense>

union ReadData {
    uint64_t key;
    uint16_t reads[4];

    ReadData() {
    }

    ReadData(uint64_t k) {
        key = k;
    }

    ReadData(ReadData &&other) : key(other.key) {
//        std::cout << "ReadData move constructor" << std::endl;
    }

    ReadData &operator=(ReadData &&other) {
//        std::cout << "ReadData move assignment" << std::endl;
        key = other.key;
        return *this;
    }

    ReadData(const ReadData &other) = default;

    ReadData &operator=(const ReadData &other) = default;
};

typedef std::vector<ReadData> ReadDataVector;

//#pragma pack(push)  /* push current alignment to stack */
//#pragma pack(1)     /* set alignment to 1 byte boundary */
struct ModelInput{// Can probably stand to lose this, started out more complex..

    ReadDataVector all_reads;
    uint32_t site_index;
    uint16_t reference;


    ModelInput() :  site_index(-1), reference(-1){
    }

    ModelInput(const ModelInput& other) = default;
    ModelInput& operator=(const ModelInput& other) = default;

    ModelInput(ModelInput&& other) :site_index(other.site_index), reference(other.reference) {
//        std::cout << "ModelInput move constructor" << std::endl;
        other.reference = 0;
        other.site_index = 0;
        all_reads = std::move(other.all_reads);
        all_reads.shrink_to_fit();
    }

    ModelInput& operator=(ModelInput&& other){
//        std::cout << "ModelInput move assignment" << std::endl;
        reference = other.reference;
        site_index = other.site_index;
        other.reference = 0;
        other.site_index = 0;
        all_reads = std::move(other.all_reads);
        return *this;
    }

    ModelInput(uint32_t site_index, int read_data_count) : site_index(site_index), reference(-1)  {
        all_reads.reserve(read_data_count);
        for (int i = 0; i < read_data_count; ++i) {
            all_reads.emplace_back(ReadData{0});
        }
//        all_reads = ReadDataVector(read_data_count, ReadData{0}     );
    }

//    ModelInput(uint16_t reference0, ReadDataVector &all_reads0) : ModelInput(-1, reference0, all_reads0) {
//    }

    ModelInput(uint32_t site_index, uint16_t reference0, ReadDataVector &all_reads0) : all_reads(all_reads0),
                                                                                       site_index(site_index),
                                                                                       reference(reference0) {
    }

    ~ModelInput() { };
};
//#pragma pack(pop)   /* restore original alignment from stack */


void PrintReads(ReadData data);
void PrintModelInput(ModelInput model_input);


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


int parseLine(char* line);

void printMemoryUsage(char const *string1 = "");
int getMemoryUsageVmSize();

int getMemoryUsageVmRSS();


#endif
