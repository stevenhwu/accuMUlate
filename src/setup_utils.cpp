//
// Created by steven on 5/13/15.
//

#include <algorithm/em_algorithm_thread_mutation.h>
#include <vector>
#include "setup_utils.h"

void CustomFilterGenomeData(GenomeData &genome_data) {

//    Custom filter, remove sample 47 51 44   ->  7 8 10

    int read_lower_bound = 6;
    int read_upper_bound = 150;
    for (int g = genome_data.size()-1 ; g >= 0; --g) {
        ReadDataVector &read_vector = genome_data[g].all_reads;
        RemoveDescendantRDV(read_vector, 10);
        RemoveDescendantRDV(read_vector, 8);
        RemoveDescendantRDV(read_vector, 7);

//        std::cout << "\t" << data.all_reads.size() << "\t" << data.all_reads.capacity() << "\t" << std::endl;

        int sum = 0;
        for (int j = 0; j < 4; ++j) {
            sum += read_vector[0].reads[j];
        }

        if(sum < read_lower_bound || sum > read_upper_bound){//Bad ref, remove site
//            PrintModelInput(data);
            RemoveSiteGenomeData(genome_data, g);

        }

        else {//Check each descendant
            for (int i = read_vector.size()-1 ; i > 0; --i) {
                int sum = 0;
                for (int j = 0; j < 4; ++j) {
                    sum += read_vector[i].reads[j];
                }
                if (sum < read_lower_bound || sum > read_upper_bound) {
                    RemoveDescendantRDV(read_vector, i);

//                    std::cout << g << "\t" << sum << "\t" << i << "\t" << data.all_reads.size() << "\t" << data.all_reads[i].key <<"\t" << data.all_reads[i-1].key <<"\t" << data.all_reads[i+1].key <<"\t";
                }
            }

            if (read_vector.size() ==1) {
//                std::cout << g << "\t" << data.all_reads.size() << std::endl;
//                for (int i = 0; i < 9; ++i) {
//                    PrintReads(data.all_reads[i]);
//                }

//                std::swap(genome_data[g], genome_data.back());
//                genome_data.pop_back();
                RemoveSiteGenomeData(genome_data, g);

            }
        }

        read_vector.shrink_to_fit();

    }

}




GenomeData GetGenomeData(boost::program_options::variables_map variables_map) {


    std::string file_name = variables_map["bam"].as<std::string>();// "zz_test_GenomeData_binary_subset";
//    PileupUtils::WriteGenomeDataToBinary(file_name, genome_data);
    GenomeData genome_data;
    PileupUtils::ReadGenomeDataFromBinary(file_name, genome_data);
    return genome_data;

}


uint16_t MaxElementIndex(ReadData &read_data){
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

    GenomeData genome_data = GetGenomeData(variables_map);
    printMemoryUsage("Read genomeData");

    int index[12];
    int ref_diff = 0;
    int des_diff = 0;
    std::map<int, int> counter;
    for (auto &site : genome_data) {
        bool new_line = false;
        uint16_t ref = site.reference;
        ReadDataVector &rd = site.all_reads;
        index[0] = MaxElementIndex(rd[0]);

//        if(ref != index[0]){
//            cout << "Ref: " << ref << "\t" << index[0] << "===";
//            ref_diff ++;
//            new_line = true;
//        }
        int local_count = 0;
        for (int i = 1; i < 12; ++i) {
            index[i] = MaxElementIndex(rd[i]);
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





void SummariseRealData(GenomeData &genome_data ) {



    clock_t t1;
    t1 = clock();


    int index[12];
    int ref_diff = 0;
    int des_diff = 0;
    std::map<int, int> counter;
    for (auto &site : genome_data) {
        bool new_line = false;
        uint16_t ref = site.reference;
        ReadDataVector &rd = site.all_reads;
        index[0] = MaxElementIndex(rd[0]);

//        if(ref != index[0]){
//            cout << "Ref: " << ref << "\t" << index[0] << "===";
//            ref_diff ++;
//            new_line = true;
//        }
        int local_count = 0;
        for (int i = 1; i < rd.size(); ++i) {
            index[i] = MaxElementIndex(rd[i]);
            if(index[i] != index[0] && (rd[i].reads[index[i]]>0) ){
//                cout << index[i] ;
//                cout << index[0] << index[1] << endl;
                des_diff ++;
                local_count++;
//                break;
                new_line = true;
            }
        }
        if(local_count>0){
//            ref_diff++;
//            std::cout << local_count << ":   ";
//            for (int i = 0; i < rd.size(); ++i) {
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

void RemoveDescendantRDV(ReadDataVector &read_vector, int d) {
    std::swap(read_vector[d], read_vector.back());
    read_vector.pop_back();
}

void RemoveSiteGenomeData(GenomeData &genome_data, int g) {
//    std::cout << "\tRemove g:"<<g << "\t" << genome_data[g].site_index << "\t" ;
    std::swap(genome_data[g], genome_data.back());
//    std::cout << genome_data[g].site_index  << "\t" ;
    genome_data.pop_back();
//    std::cout << genome_data[g].site_index  << std::endl;
}



