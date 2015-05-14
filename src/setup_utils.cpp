//
// Created by steven on 5/13/15.
//

#include "setup_utils.h"


void PostFilterGenomeData(GenomeData &genome_data) {

//    Custom filter, remove sample 47 51 44   ->  7 8 10
//  Alter the size of each descendant
    for (int g = genome_data.size()-1 ; g >= 0; --g) {

        ReadDataVector &data_vector = genome_data[g].all_reads;
//        PrintModelInput(data);
//        std::cout << "\t" << data_vector.size() << "\t" << data_vector.capacity() << "\t";
//        data_vector.erase(data_vector.begin()+10);
//        data_vector.erase(data_vector.begin()+7, data_vector.begin()+9);
//        data_vector.shrink_to_fit();


        std::swap(data_vector[10], data_vector.back());
        data_vector.pop_back();
        std::swap(data_vector[8], data_vector.back());
        data_vector.pop_back();
        std::swap(data_vector[7], data_vector.back());
        data_vector.pop_back();
//        std::cout << "\t" << data_vector.size() << "\t" << data_vector.capacity() << "\t" << std::endl;

/*
        int sum = 0;
        for (int j = 0; j < 4; ++j) {
            sum += data_vector[0].reads[j];
        }

        if(sum < 6){//Bad ref, remove site
//            PrintModelInput(data);
//            genome_data.erase(genome_data.begin()+g);
            std::swap(genome_data[g], genome_data.back());
            genome_data.pop_back();
        }
//        else if(sum>150){
//            PrintReads(data_vector[0]);
//        }
        else {//Check each descendant
//            std::cout <<g << "\t" << data_vector.capacity() << "\t" << data_vector.size() << std::endl;
            for (int i = data_vector.size()-1 ; i > 0; --i) {
                int sum = 0;
                for (int j = 0; j < 4; ++j) {
                    sum += data_vector[i].reads[j];
                }
                if (sum < 6) {
//                    std::cout << g << "\t" << sum << "\t" << i << "\t" << data_vector.size() << "\t" << data_vector[i].key <<"\t" << data_vector[i-1].key <<"\t" << data_vector[i+1].key <<"\t";
//                    data_vector.erase(data_vector.begin()+i);//FIXME: current model assume equal des, fail if do this
                    std::swap(data_vector[i], data_vector.back());
                    data_vector.pop_back();

//                    if (i != pList.size() - 1)
//                    {
//                        // Beware of move assignment to self
//                        // see http://stackoverflow.com/questions/13127455/
//                        pList[i] = std::move(pList.back());
//                    }
//                    pList.pop_back();
                }
//                else if(sum>200){
//                    PrintReads(data_vector[i]);
//                }
            }

            if (data_vector.size() ==1) {
//                std::cout << g << "\t" << data_vector.size() << std::endl;
//                for (int i = 0; i < 9; ++i) {
//                    PrintReads(data_vector[i]);
//                }

                std::swap(genome_data[g], genome_data.back());
                genome_data.pop_back();

            }
        }
*/
    }

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
        index[0] = max_element_index(rd[0]);

//        if(ref != index[0]){
//            cout << "Ref: " << ref << "\t" << index[0] << "===";
//            ref_diff ++;
//            new_line = true;
//        }
        int local_count = 0;
        for (int i = 1; i < rd.size(); ++i) {
            index[i] = max_element_index(rd[i]);
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

