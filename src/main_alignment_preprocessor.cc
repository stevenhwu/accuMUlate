
#include <iostream>
#include <map>
#include <vector>

#include <google/profiler.h>

#include <boost/program_options.hpp>
#include <time.h>
#include <algorithm/em_model_binomial.h>
#include <algorithm/em_data_binomial.h>
#include <algorithm/em_algorithm_binomial.h>
#include <stdint.h>
//#include <sys/socket.h>


#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"
#include "parsers.h"
#include "sequence_prob.h"
#include "VariantVisitor.h"
#include "variant_visitor_two.h"
#include "evolution_models/JC69.h"
#include "evolution_models/F81.h"
#include "site_prob.h"
#include "algorithm/em_algorithm_mutation.h"
#include "algorithm/em_data_mutation.h"
#include "genome_data_stream.h"

#include "boost_input_utils.h"
#include "pileup_utils.h"

using namespace std;
using namespace BamTools;


void SummariseReadsData(GenomeData base_counts) {
    size_t site_count = base_counts.size();
    for (size_t i = 0; i < site_count; ++i) {

        int sum = 0;
        int total = 0;
        for (size_t j = 0; j < base_counts[i].all_reads.size(); ++j) {

//            base_counts[i].all_reads[3]
            ReadData &reference = base_counts[i].all_reads[j];
            auto reads = reference.reads;
            uint16_t max = *std::max_element(reads, reads + 4);
            uint16_t ref_base_count = reads[base_counts[i].reference];
            uint16_t diff = abs(max - ref_base_count);
            sum += diff;
            total += ref_base_count;
//            cout << diff << " -- " << max << " " << ref_base_count << "\t== ";
//            SequenceProb::printReadData(base_counts[i].all_reads[j]) ;

        }
        double prop = (double) sum / total;
        if (prop > 0.1) {
            cout << "======= Site: " << i << " R: " << base_counts[i].reference;// << endl;

            cout << "\t==" << sum << " " << total << " " << prop << endl;

            for (size_t j = 0; j < base_counts[i].all_reads.size(); ++j) {
//                SequenceProb::printReadData(base_counts[i].all_reads[j]);
            }
        }
//        cout << "================================="<< endl;
    }
//    cout << "================Done: SequenceProb. Total: " << site_count << endl;


    /*
======= Site: 64 R: 1	==55 23 2.3913
======= Site: 163 R: 3	==134 15 8.93333
======= Site: 805 R: 3	==25 47 0.531915
======= Site: 808 R: 2	==18 49 0.367347
======= Site: 884 R: 1	==171 28 6.10714
======= Site: 969 R: 1	==255 2 127.5
======= Site: 4698 R: 3	==12 94 0.12766
======= Site: 8622 R: 0	==45 98 0.459184
======= Site: 8626 R: 2	==25 103 0.242718
======= Site: 9415 R: 3	==232 21 11.0476
======= Site: 9436 R: 0	==188 45 4.17778
======= Site: 9459 R: 2	==156 80 1.95
======= Site: 9473 R: 3	==114 117 0.974359

    */

}

int main(int argc, char** argv){

    boost::program_options::variables_map variable_map;

    BoostUtils::ParseCommandLinkeInput(argc, argv, variable_map);

    GenomeData genome_data;
//    PileupUtils::CreatePileupAlignment(variable_map, genome_data, 1);
    PileupUtils::CreatePileupAlignment(variable_map, genome_data, 2);

    clock_t start0 = clock();
    std::string file_name = "zz_test_GenomeData_binary_subset";
    PileupUtils::WriteGenomeDataToBinary(file_name, genome_data);

//    PileupUtils::ReadGenomeDataFromBinary(file_name, genome_data);

    printf("Total ReadWrite Time: %f \n", ((double) (clock()-start0)/ CLOCKS_PER_SEC) );
//    SummariseReadsData(genome_datas);



}




