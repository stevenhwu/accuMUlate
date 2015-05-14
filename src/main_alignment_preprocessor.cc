
#include <iostream>
#include <io_data/boost_input_utils.h>
#include "io_data/pileup_utils.h"



/*
cd build
cmake ..
make AlignmentPreprocessor


cd data
../build/src/AlignmentPreprocessor   -b tt_subset.bam -r tt-ref.fasta -c params.ini  --output_binary_file test_interval_1 -i interval_1.bed
../build/src/AlignmentPreprocessor   -b tt_subset.bam -r tt-ref.fasta -c params.ini  --output_binary_file test_all

python ../src/utils/split_beds.py tt_subset.bam 3 .
*/


int main(int argc, char** argv){

    boost::program_options::variables_map variable_map;

    BoostUtils::ParseCommandLinkeInput(argc, argv, variable_map);

    GenomeData genome_data;
//    PileupUtils::CreatePileupAlignment(variable_map, genome_data, 1);
    PileupUtils::CreatePileupAlignment(variable_map, genome_data, 2);

//    clock_t start0 = clock();
//    std::string file_name = "test_GenomeData_binary_subset";
//    PileupUtils::WriteGenomeDataToBinary(file_name, genome_data);
//    PileupUtils::ReadGenomeDataFromBinary(file_name, genome_data);

//    printf("Total ReadWrite Time: %f \n", ((double) (clock()-start0)/ CLOCKS_PER_SEC) );
//    SummariseReadsData(genome_datas);



}




