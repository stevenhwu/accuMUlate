//
// Created by steven on 5/13/15.
//

#include "setup_utils.h"


void PostFilterGenomeData(GenomeData &genome_data) {

//    Custom filter, remove sample 47 51 44   ->  7 8 10
//
//    for (auto &data : genome_data) {
    for (int g = genome_data.size()-1 ; g >= 0; --g) {

        ModelInput &data = genome_data[g];
//        for (int i = 0; i < data.all_reads.size(); ++i) {
//            std::cout << i << ":" << data.all_reads[i].key << " ";
//        }
//        std::cout << "" << std::endl;
        data.all_reads.erase(data.all_reads.begin()+10);
        data.all_reads.erase(data.all_reads.begin()+7, data.all_reads.begin()+9);

//        for (int i = 0; i < data.all_reads.size(); ++i) {
//            std::cout << i << ":" << data.all_reads[i].key << " ";
//        }
//        std::cout << "" << std::endl;
//        data.all_reads.erase(data.all_reads.begin()+11);
//        data.all_reads.erase(data.all_reads.begin()+9);
//        data.all_reads.erase(data.all_reads.begin()+1, data.all_reads.begin()+7);

        int sum = 0;

        for (int j = 0; j < 4; ++j) {
            sum += data.all_reads[0].reads[j];
        }
        if(sum < 6){
//            data.all_reads.clear();
//            std::cout << g << ":" << genome_data.size() << " ";
//
//            std::cout << genome_data.size() << "\t" << data.all_reads[0].key << "\t" << sum << std::endl;
//            for (int j = 0; j < 4; ++j) {
//                std::cout << data.all_reads[0].reads[j]  << std::endl;
//            }

            genome_data.erase(genome_data.begin()+g);
//            break;

        }
        else {
            int zero_count = 0;
            for (int i = data.all_reads.size() - 1; i > 0; --i) {
                int sum = 0;

                for (int j = 0; j < 4; ++j) {
                    sum += data.all_reads[i].reads[j];
                }
//            std::cout << i << "\t" << data.all_reads[i].key << "\t" << sum <<std::endl;
//            std::cout << i << ":" << sum << " ";
                if (sum < 6) {
//                data.all_reads.erase(data.all_reads.begin()+i);//FIXME: current model assume equal des, fail if do this
                    data.all_reads[i] = 0;
                    zero_count ++;
                }

            }
//        std::cout << "" << std::endl;
//        for (int i = 0; i < data.all_reads.size(); ++i) {
//            std::cout << i << ":" << data.all_reads[i].key << " ";
//        }
//        std::cout << "" << std::endl;
//        std::cout << "FINAL SIZE: " << data.all_reads.size() << std::endl;
//        std::exit(13);
            if (zero_count == data.all_reads.size() -1) {
//                data.all_reads.clear();
//                std::cout << "FINAL SIZE: " << data.all_reads.size() << "\t" << zero_count << std::endl;
//                for (int i = 0; i < data.all_reads.size(); ++i) {
//                    std::cout << i << ":" << data.all_reads[i].key << " ";
//                }
//                 std::cout << "" << std::endl;

                genome_data.erase(genome_data.begin() + g);

            }
        }
//        for (int i = 0; i < data.all_reads.size(); ++i) {
//            int sum = 0;
//            for (int j = 0; j < 4; ++j) {
//                sum += data.all_reads[i].reads[j];
//            }
////            if(sum<6){
////                std::cout << i << "\t" ;
////                for (int j = 0; j < 4; ++j) {
////                    std::cout << data.all_reads[i].reads[j] << " "  ;
////                }
////                std::cout << "" << std::endl;
////            }
////            std::cout << i << "\t" << data.all_reads[i].key << std::endl;
//        }
//        std::cout << "" << std::endl;
//        std::cout << "\n=======\n" << std::endl;
//        exit(32);
    }
//FinalSummary: 5.57299e-01	6.02335e-06	1.32571e-04	9.99867e-01	-238817.82286	7.70999e-11
//EM Clock time: 32	32148156

}


GenomeData getGenomeData(boost::program_options::variables_map variables_map) {


    std::string file_name = variables_map["bam"].as<std::string>();// "zz_test_GenomeData_binary_subset";
//    PileupUtils::WriteGenomeDataToBinary(file_name, genome_data);
    GenomeData genome_data;
    PileupUtils::ReadGenomeDataFromBinary(file_name, genome_data);
    return genome_data;

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
        if(local_count>3){
            ref_diff++;
//            for (int i = 0; i < 12; ++i) {
//                std::cout << index[i] << ":" << rd[i].reads[index[i]] << " ";
//            }
//            std::cout << "" << std::endl;
        }
        counter[local_count]++;
//        if(new_line) {
//            cout << "\t" <<  local_count << "\t" << endl;
//        }

    }
    std::cout << "Summary: " << ref_diff << "\t" << des_diff << "\t" << genome_data.size() << std::endl;
    for (auto item : counter) {
        std::cout << item.first << "\t" << item.second << std::endl;
    }

}


