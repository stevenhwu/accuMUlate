
#include <iostream>
#include <map>
#include <vector>
#include <boost/program_options.hpp>
#include <random>
#include <algorithm/em_algorithm_thread_mutation.h>
#include <algorithm/em_algorithm_mutation.h>

#include "api/BamReader.h"
#include "evolution_models/F81.h"
#include "boost_input_utils.h"
#include "pileup_utils.h"
#include "algorithm/em_algorithm_mutation.h"
#include "algorithm/em_algorithm_mutation_v1.h"
#include "sequencing_factory.h"

#include "sys/types.h"
#include "sys/sysinfo.h"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
//#include "em_algorithm.h"
#include "mutation_model_multi_categories.h"

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wdeprecated"


using namespace std;
using namespace BamTools;
/*
 *
 * Custom filter, remove sample 47 51 44   ->  7 8 10
M50_CAGATC_L001	9
M47_ATTGAGGA_L007	8
M44_GCCAAT_L001	7
M531_CTTGTA_L005	11
M50_CAGATC_L005	9
M40_ACAGTG_L005	6
M44_GCCAAT_L005	7
M40_ACAGTG_L001	6
M531_CTTGTA_L001	11
M29_CAGATCTG_L007	5
M28_TGACCA_L005	4
M25_GTGTTCTA_L007	3
M28_TGACCA_L001	4
M19_AACGTGAT_L007	1
M0_CGATGT_L005	0
M51_GGAGAACA_L007	10
M20_AAACATCG_L007	2
M0_CGATGT_L001	0
 */

//namespace po = boost::program_options;
int RunBasicProbCalc(GenomeData base_counts, ModelParams params);

void testCalWeighting(MutationProb &mutation_prob, std::vector<SequenceProb> sp);

void AddSimulatedData(ModelParams &params, std::vector<SequenceProb> &sp, int descendant_count, size_t fake_sample_count, double fake_prop);

void SummariseReadsData(GenomeData base_counts);


void CreateMutationModel(MutationModel &model, GenomeData &genome_data, ModelParams &params);

GenomeData getGenomeData(boost::program_options::variables_map variables_map);

MutationModelMultiCategories CreateMutationModelMulti(GenomeData genome_data, ModelParams params);

void RunEmWithRealDataMultiThread(boost::program_options::variables_map map, size_t thread_count);

void SummariseRealData(boost::program_options::variables_map map);

void PostFilterGenomeData(GenomeData &genome_data);

void RunEmWithRealDataOriginal(boost::program_options::variables_map variables_map) {
    ModelParams params = BoostUtils::CreateModelParams(variables_map);
    MutationProb mutation_prob = MutationProb(params);
    F81 evo_model0(mutation_prob);
    GenomeData genome_data = getGenomeData(variables_map);

    cout << "\n========================\nStart em_algorithm:" << endl;
    clock_t t_start, t_end;
    t_start = clock();

    std::vector<SequenceProb> sp;
    SequenceProb::CreateSequenceProbV1(sp, genome_data, params);
    std::vector<std::unique_ptr<EmData>> em_site_data;
    for (auto &seq_prob: sp) {
        em_site_data.emplace_back(new EmDataMutationV1(seq_prob, evo_model0));
    }
    EmModelMutationV1 em_model0(evo_model0);
    t_end = clock();
    cout << "Time old preprocess Data: " << (t_end - t_start) / CLOCKS_PER_SEC << "\t" << (t_end - t_start) << endl;

    t_start = clock();
    EmAlgorithmMutationV1 em_alg2(2, em_site_data, em_model0);
    em_alg2.RunEM();
    t_end = clock();
    cout << "Time old: " << (t_end - t_start) / CLOCKS_PER_SEC << "\t" << (t_end - t_start) << endl;

    em_alg2.PrintSummary();

}
//void RunEmWithRealData(GenomeData &genome_data, ModelParams params) {
void RunEmWithRealDataOneStep(boost::program_options::variables_map variables_map) {

    ModelParams params = BoostUtils::CreateModelParams(variables_map);
    MutationProb mutation_prob = MutationProb(params);
    F81 evo_model0(mutation_prob);
    MutationModel mutation_model = MutationModel(evo_model0);

    printMemoryUsage("Start ");
    clock_t t1;
    t1 = clock();


    std::string file_name = variables_map["bam"].as<string>();// "zz_test_GenomeData_binary_subset";
    GenomeDataStream gd_stream_read = GenomeDataStream( file_name, false);

    uint64_t total_base_count2 = 0;
    uint64_t sequence_count2 = 0;
    gd_stream_read.ReadHeader(total_base_count2, sequence_count2);
    std::cout << "Sequence count: " << sequence_count2 << " Bases: " << total_base_count2 << std::endl;

//    gd_stream_read.ReadGenomeData();
//    GenomeData genome_data;
//    getGenomeData(genome_data, variables_map);
    printMemoryUsage("Read genomeData");
//    cout << "===== Init site_count: " << genome_data.size() << endl;
//    cout << "===== Setup EmData:" << endl;
//
//    cout << "Time: read genome data: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;
    t1 = clock();


    std::vector<SiteGenotypesIndex> sgi;
    sgi.reserve(total_base_count2);
    printMemoryUsage("Before factory");
    SequencingFactory sequencing_factory(params);
    printMemoryUsage("init factory");


    for (int i = 0; i < total_base_count2; ++i) {
        ModelInput m = gd_stream_read.ReadOne();
        sequencing_factory.CreateSequenceProbsVector(sgi, m);
    }

    cout << "Time init seq latest: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;
    printMemoryUsage("after factory");

    gd_stream_read.close(); //close gd_stream


//    MutationModel mutation_model = MutationModel(evo_model0);

    MutationModel::AddGenotypeFactory(sequencing_factory);
    printMemoryUsage("add genotypes");
    std::cout << sgi.size() << std::endl;
    mutation_model.MoveSequenceProb(std::move(sgi));
    printMemoryUsage("add seq probs");
    std::cout << sgi.size() << std::endl;

    cout << "===== Done preprocess. Final site count: " << mutation_model.GetSiteCount() << endl;
    printMemoryUsage("Created Mutation Model");

    cout << "===== Setup EM" << endl;
    std::vector<std::unique_ptr<EmModel>> em_model2;
    em_model2.emplace_back(new EmModelMutation(mutation_model));
    em_model2.emplace_back(new EmModelMutation(mutation_model));
    printMemoryUsage("Double object?");

    cout << "\n========================\nStart em_algorithm:" << endl;
    clock_t t_start, t_end;
    t_start = clock();
    EmAlgorithmMutation em_alg0(em_model2);

    printMemoryUsage("EM");
//        exit(2);

    em_alg0.RunEM();
    em_alg0.PrintSummary();
    t_end = clock();
    cout << "Time new: " << (t_end - t_start) / CLOCKS_PER_SEC << "\t" << (t_end - t_start) << endl << endl;
    printMemoryUsage();


//    exit(2);
//    t_start = clock();
//    std::vector<SequenceProb> sp;
//    SequenceProb::CreateSequenceProbV1(sp, genome_data, params);
//    std::vector<std::unique_ptr<EmData>> em_site_data;
//    for (auto &seq_prob: sp) {
//        em_site_data.emplace_back(new EmDataMutationV1(seq_prob, evo_model0));
//    }
//    EmModelMutationV1 em_model0(evo_model0);
//    t_end = clock();
//    cout << "Time old preprocess Data: " << (t_end - t_start) / CLOCKS_PER_SEC << "\t" << (t_end - t_start) << endl;
//
//    t_start = clock();
//    EmAlgorithmMutationV1 em_alg2(2, em_site_data, em_model0);
//    em_alg2.Run();
//    t_end = clock();
//    cout << "Time old: " << (t_end - t_start) / CLOCKS_PER_SEC << "\t" << (t_end - t_start) << endl;
//
//    em_alg2.PrintSummary();
//    em_alg0.PrintSummary();
}

void RunEmWithRealData(boost::program_options::variables_map variables_map) {

    printMemoryUsage("Start ");
    ModelParams params = BoostUtils::CreateModelParams(variables_map);
    MutationProb mutation_prob (params);
    F81 evo_model0(mutation_prob);
    MutationModel mutation_model (evo_model0);

    printMemoryUsage("Start Basic");

    clock_t t1;
    t1 = clock();

    GenomeData genome_data = getGenomeData(variables_map);
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

    GenomeData genome_data = getGenomeData(variables_map);
    PostFilterGenomeData(genome_data);
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

void PostFilterGenomeData(GenomeData &genome_data) {

//    Custom filter, remove sample 47 51 44   ->  7 8 10
//
    for (auto &data : genome_data) {
        for (int i = 0; i < data.all_reads.size(); ++i) {
//            std::cout << i << "\t" << data.all_reads[i].key << std::endl;
        }
        data.all_reads.erase(data.all_reads.begin()+10);
        data.all_reads.erase(data.all_reads.begin()+7, data.all_reads.begin()+9);

//     std::cout << data.all_reads.size() << std::endl;
        for (int i = 0; i < data.all_reads.size(); ++i) {
//            std::cout << i << "\t" << data.all_reads[i].key << std::endl;
        }
//        std::cout << "\n=======\n" << std::endl;
//        exit(32);
    }
//FinalSummary: 5.57299e-01	6.02335e-06	1.32571e-04	9.99867e-01	-238817.82286	7.70999e-11
//EM Clock time: 32	32148156

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


void AddSimulatedData(ModelParams &params, std::vector<SequenceProb> &sp, int descendant_count, size_t fake_sample_count, double fake_prop) {

    random_device rd;
    mt19937 e2(rd());
    uniform_int_distribution<uint16_t> uniform_dist(0, 5);
    uniform_int_distribution<uint16_t> uniform3(0, 3);

    size_t fake_diff_count = fake_sample_count * fake_prop;
    cout << "========= Add simulated data:" << fake_sample_count << " with fake_diff_count: " << fake_diff_count << endl;

    descendant_count++;

    for (size_t s = 0; s < fake_sample_count; ++s) {

        ReadDataVector bcalls(descendant_count, ReadData{0});
        for (int i = 0; i < descendant_count; ++i) {
            bcalls[i].key=0;
            for (int j = 0; j < 4; ++j) {
                bcalls[i].reads[j] = uniform_dist(e2);
            }
            bcalls[i].reads[0] = 20 + uniform_dist(e2);
        }

        if(s < fake_diff_count){ //diff
            bcalls[0].reads[0] = uniform_dist(e2);
            bcalls[0].reads[3] = 20 + uniform_dist(e2);
            for (int i = 1; i < descendant_count; ++i) {
                bcalls[i].reads[uniform3(e2)] = 20 + uniform_dist(e2);
            }
        }
        uint16_t ref_index = 0;
        sp.emplace_back(SequenceProb( ModelInput{ref_index, bcalls} , params) );
    }

}



GenomeData getGenomeData(boost::program_options::variables_map variables_map) {


    std::string file_name = variables_map["bam"].as<string>();// "zz_test_GenomeData_binary_subset";
//    PileupUtils::WriteGenomeDataToBinary(file_name, genome_data);
    GenomeData genome_data;
    PileupUtils::ReadGenomeDataFromBinary(file_name, genome_data);
    return genome_data;

}

int main(int argc, char** argv){
    Eigen::initParallel();//TODO: triple check this!!

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    boost::program_options::variables_map variables_map;
    BoostUtils::ParseCommandLinkeInput(argc, argv, variables_map);

    SummariseRealData(variables_map);
//    exit(83);
    {
//        RunEmWithRealData(genome_data, params);
//        RunEmWithRealDataOneStep(variables_map);

//        RunEmWithRealData(variables_map);
        int thread_count = variables_map["thread"].as<int>();
        RunEmWithRealDataMultiThread(variables_map, thread_count);

//        RunEmWithRealDataOriginal(variables_map);

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

uint16_t max_element_index(ReadData &read_data){
    uint16_t max = *(  std::max_element(read_data.reads, read_data.reads+4) );
    for (int i = 0; i < 4; ++i) {
        if(max==read_data.reads[i])
            return i;
    }
    return -1;
}
void SummariseRealData(boost::program_options::variables_map variables_map) {


    printMemoryUsage("Start ");
    ModelParams params = BoostUtils::CreateModelParams(variables_map);
//    MutationProb mutation_prob (params);
//    F81 evo_model0(mutation_prob);
//    MutationModel mutation_model (evo_model0);

    printMemoryUsage("Start Basic");

    clock_t t1;
    t1 = clock();

    GenomeData genome_data = getGenomeData(variables_map);
    printMemoryUsage("Read genomeData");

    int index[12];
    int ref_diff = 0;
    int des_diff = 0;
    std::map<int, int> counter;
    for (auto &site : genome_data) {
        bool new_line = false;
        uint16_t ref = site.reference;
        ReadDataVector &rd = site.all_reads;
        index[0] = max_element_index(rd[0]);

//        if(ref != index[0]){
//            cout << "Ref: " << ref << "\t" << index[0] << "===";
//            ref_diff ++;
//            new_line = true;
//        }
        int local_count = 0;
        for (int i = 1; i < 12; ++i) {
            index[i] = max_element_index(rd[i]);
            if(index[i] != index[0]){
//                cout << index[i] ;
//                cout << index[0] << index[1] << endl;
                des_diff ++;
                local_count++;
//                break;
                new_line = true;
            }
        }
        if(local_count>4){
            ref_diff++;
        }
        counter[local_count]++;
//        if(new_line) {
//            cout << "\t" <<  local_count << "\t" << endl;
//        }

    }
    cout << "Summary: " << ref_diff << "\t" << des_diff << "\t" << genome_data.size() << endl;
    for (auto item : counter) {
        cout << item.first << "\t" << item.second << endl;
    }

}


