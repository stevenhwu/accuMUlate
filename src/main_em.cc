

#include <iostream>
#include <map>
#include <vector>
#include <boost/program_options.hpp>
#include <algorithm/em_algorithm_thread_mutation.h>
#include <algorithm/em_algorithm_mutation.h>
#include <io_data/boost_input_utils.h>
#include <io_data/pileup_utils.h>
#include <evolution_models/F81.h>
#include "setup_utils.h"


#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wdeprecated"


using namespace std;
using namespace BamTools;
/*
 *
 * Custom filter, remove sample 47 51 44   ->  7 8 10
M47_ATTGAGGA_L007	8
M44_GCCAAT_L001	7
M51_GGAGAACA_L007	10

 M50_CAGATC_L001	9
 M50_CAGATC_L005	9
 M531_CTTGTA_L001	11
 M531_CTTGTA_L005	11
 M40_ACAGTG_L005	6
 M40_ACAGTG_L001	6
 M44_GCCAAT_L005	7

 M29_CAGATCTG_L007	5
 M28_TGACCA_L001	4
 M28_TGACCA_L005	4

 M25_GTGTTCTA_L007	3
 M20_AAACATCG_L007	2
 M19_AACGTGAT_L007	1
 M0_CGATGT_L001	0
 M0_CGATGT_L005	0



 */



void CreateMutationModel(MutationModel &model, GenomeData &genome_data, ModelParams &params);

MutationModelMultiCategories CreateMutationModelMulti(GenomeData genome_data, ModelParams params);

void RunEmWithRealDataMultiThread(boost::program_options::variables_map map, size_t thread_count);


void RunEmWithRealData(boost::program_options::variables_map variables_map) {

    printMemoryUsage("Start ");
    ModelParams params = BoostUtils::CreateModelParams(variables_map);
    MutationProb mutation_prob (params);
    F81 evo_model0(mutation_prob);
    MutationModel mutation_model (evo_model0);

    printMemoryUsage("Start Basic");

    clock_t t1;
    t1 = clock();

    GenomeData genome_data = GetGenomeData(variables_map);
    CustomFilterGenomeData(genome_data);//TODO: Add this back

    printMemoryUsage("Read genomeData");

    cout << "Time: read genome data: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;
    cout << "===== Setup EmData. Init site_count: " << genome_data.size() << endl;

    t1 = clock();
    CreateMutationModel(mutation_model, genome_data, params);
    cout << "Time Create Mutation Model: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;

    cout << "===== Done preprocess. Final site count: " << mutation_model.GetSiteCount() << endl;
    printMemoryUsage("Created Mutation Model");
    cout << "===== Setup EM" << endl;

    std::vector<std::unique_ptr<EmModel>> em_model2;
    em_model2.emplace_back(new EmModelMutation(mutation_model));
    em_model2.emplace_back(new EmModelMutation(mutation_model));
    const string &outfile_prefix = variables_map["outfile"].as<string>();

    clock_t t_start, t_end;
    t_start = clock();

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    EmAlgorithmMutation em_alg0(em_model2);
    em_alg0.SetOutfilePrefix(outfile_prefix);
    printMemoryUsage("EM");
    em_alg0.Run();
    em_alg0.PrintSummary();
    std::string summary = em_alg0.GetEMSummary();
    std::cout << summary << std::endl;

    t_end = clock();
    cout << "EM Clock time: " << (t_end - t_start) / CLOCKS_PER_SEC << "\t" << (t_end - t_start) << endl << endl;

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "EM started at "<< std::ctime(&start_time) << "EM finished at " << std::ctime(&end_time)  << "EM elapsed time: " << elapsed_seconds.count() << "s\n";

    printMemoryUsage("End EM");

}


void RunEmWithRealDataMultiThread(boost::program_options::variables_map variables_map, size_t thread_count) {

    printMemoryUsage("Start ");
    ModelParams params = BoostUtils::CreateModelParams(variables_map);
    MutationProb mutation_prob (params);
    F81 evo_model0(mutation_prob);
//    MutationModel mutation_model (evo_model0);

    printMemoryUsage("Start Basic");

    clock_t t1;
    t1 = clock();

    GenomeData genome_data = GetGenomeData(variables_map);
    SummariseRealData(genome_data);
    CustomFilterGenomeData(genome_data);//TODO: Add this back
    SummariseRealData(genome_data);
    printMemoryUsage("Read genomeData");

    cout << "Time: read genome data: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;
    cout << "===== Setup EmData. Init site_count: " << genome_data.size() << endl;

    t1 = clock();

    SequencingFactory sequencing_factory(params);
    printMemoryUsage("init factory");
    sequencing_factory.CreateSequenceProbsVector(genome_data);

    MutationModelMultiCategories modelMulti (2, evo_model0, sequencing_factory);

    cout << "Time Create Mutation Model: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;

    cout << "===== Done preprocess. Final site count: " << modelMulti.GetSiteCount() << endl;
    printMemoryUsage("Created Mutation Model");
    cout << "======= Setup EM ======" << endl;

    const string &outfile_prefix = variables_map["outfile"].as<string>();

    clock_t t_start, t_end;
    t_start = clock();

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    EmAlgorithmThreadMutation em_algM(modelMulti, thread_count);
    em_algM.SetOutfilePrefix(outfile_prefix);
    printMemoryUsage("EM");
    em_algM.Run();
    em_algM.PrintSummary();
    std::string summary = em_algM.GetEMSummary();
    std::cout << "FinalSummary: " << summary << std::endl;


    t_end = clock();
    cout << "EM Clock time: " << (t_end - t_start) / CLOCKS_PER_SEC << "\t" << (t_end - t_start) << endl << endl;

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "EM started at "<< std::ctime(&start_time) << "EM finished at " << std::ctime(&end_time)  << "EM elapsed time: " << elapsed_seconds.count() << "s\n";

    printMemoryUsage("End EM");

}

MutationModelMultiCategories CreateMutationModelMulti(GenomeData genome_data, ModelParams params) {


    clock_t t1;
    t1 = clock();


//    std::vector<SiteGenotypesIndex> sgi;
    printMemoryUsage("Before factory");
    SequencingFactory sequencing_factory(params);
    printMemoryUsage("init factory");

    t1 = clock();
    sequencing_factory.CreateSequenceProbsVector(genome_data);
    cout << "Time init seq latest: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;
    printMemoryUsage("after factory");


    MutationProb mutation_prob (params);
    F81 evo_model0(mutation_prob);


//    MutationModel mutation_model = MutationModel(evo_model0);
    MutationModelMultiCategories model (2, evo_model0, sequencing_factory);

//    MutationModel::AddGenotypeFactory(sequencing_factory);
//    std::cout << "SiteGenotypeIndex Size: " << sgi.size() << std::endl;
//    printMemoryUsage("add genotypes");
//    mutation_model.MoveSequenceProb(std::move(sgi));
//    printMemoryUsage("add seq probs");
    return model;

}

void CreateMutationModel(MutationModel &mutation_model, GenomeData &genome_data, ModelParams &params) {

    printMemoryUsage("Before factory");

    std::vector<SiteGenotypesIndex> sgi;
    SequencingFactory sequencing_factory(params);
    sequencing_factory.CreateSequenceProbsVector(sgi, genome_data);
    printMemoryUsage("after factory");

    MutationModel::AddGenotypeFactory(sequencing_factory);
    std::cout << "SiteGenotypeIndex Size: " << sgi.size() << std::endl;
    printMemoryUsage("Added genotypes");
    mutation_model.MoveSequenceProb(std::move(sgi));
    printMemoryUsage("Added seq probs");


}

struct ZZ{// Can probably stand to lose this, started out more complex..


    ReadDataVector all_reads;
    uint32_t site_index;
    uint16_t reference;

//
};



int main(int argc, char** argv){
    Eigen::initParallel();//TODO: triple check this!!

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    boost::program_options::variables_map variables_map;
    BoostUtils::ParseCommandLinkeInput(argc, argv, variables_map);


//    SummariseRealData(variables_map);

    {
        int thread_count = variables_map["thread"].as<int>();
//        RunEmWithRealData(genome_data, params);

//        RunEmWithRealData(variables_map);
        RunEmWithRealDataMultiThread(variables_map, thread_count);


    }


    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Program started at "<< std::ctime(&start_time) <<
            "Program finished at " << std::ctime(&end_time)  <<
            "Total elapsed time: " << elapsed_seconds.count() << "s" << std::endl;


    printMemoryUsage("End");

}
