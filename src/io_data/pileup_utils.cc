/*
 * pileup_utils.cc
 *
 *  Created on: 12/13/14
 *      Author: Steven Wu
 */


#include <stdio.h>
#include "pileup_utils.h"

namespace PileupUtils{

    uint64_t print_count = 1000000;
    float fix_counter = 1e6;
    float time_scaler = fix_counter/print_count;




    void CreatePileupAlignment(boost::program_options::variables_map &vm, GenomeData &genome_data,
            int variant_visitor_index) {

        BamReader experiment;
        RefVector references;
        SamHeader header;
        LocalBamToolsUtils::Fasta reference_genome;

        BoostUtils::ExtractInputVariables(vm, genome_data, experiment, references, header, reference_genome);

        LocalBamToolsUtils::PileupEngine pileup;

        if(variant_visitor_index == 1){
            CreatePileupV1(vm, genome_data, references, header, reference_genome, pileup);
        }
        else if(variant_visitor_index == 2) {
            CreatePileupV2(vm, genome_data, references, header, reference_genome, pileup);
        }

        BamAlignment ali;

        uint64_t ali_counter = 0;
        clock_t start = clock();
        clock_t time_stored =start;

        global_count[0] = 0;global_count[1] = 0;global_count[2] = 0;

        std::cerr.setstate(std::ios_base::failbit) ; //Supress bamtool camplain, "Pileup::Run() : Data not sorted correctly!"
        //  TODO: Only allocate interval-sized memory vector
        if (vm.count("intervals")){
            BedFile bed (vm["intervals"].as<std::string>());
            BedInterval region;
            while(bed.get_interval(region) == 0){
                int ref_id = experiment.GetReferenceID(region.chr);
                experiment.SetRegion(ref_id, region.start, ref_id, region.end);

                while( experiment.GetNextAlignment(ali) ){
                    pileup.AddAlignment(ali);
                    ali_counter++;
                    VerboseAlignmentInfo(time_stored, ali_counter);
                }
            }
        }
        else{
            while( experiment.GetNextAlignment(ali)){           // Fast, 0.2s
                pileup.AddAlignment(ali);                     // AddAlignment ~2s + visitor ~3s
                ali_counter++;
                VerboseAlignmentInfo(time_stored, ali_counter);
            }
            pileup.Flush();
//            std::cout << "====DEBUG: G0: " << global_count[0] << "\tG1: " << global_count[1] << "\tG2: " << global_count[2] << std::endl;

/*      with -O0
        Let's work with 50k count for now, 4.5~5s   6788 base_count, 2s overhead
        //Init test: ~1 to 1.5s per 10k.
        // ~25s for 30.7k base_count, 194927 ali_count
        // 10s for 194927 count without visitor
        update: 50k ~ 3.3~4s
        Full 194927 ali_count: 15-18s
        Full EM: ~130s

*/
        }

        printf("Total time V%d: %fs. Ali_Count:%lu reads_count:%d %d %d genome_data.size():%lu\n", variant_visitor_index,
                ((double) (clock()-start)/ CLOCKS_PER_SEC),
                ali_counter, global_count[0], global_count[1], global_count[2], genome_data.size());

        std::cerr.clear() ;
        experiment.Close();
        reference_genome.Close();

//        streambuf *old = cout.rdbuf(); // <-- save
//        stringstream ss;
//
//        cout.rdbuf (ss.rdbuf());       // <-- redirect
//        foobar();                      // <-- call
//        cout.rdbuf (old);
    }

    void VerboseAlignmentInfo(clock_t &time_stored, uint64_t ali_counter) {
        if (ali_counter % print_count == 0){  //3800 => 67
            clock_t time_temp = clock() ;
            time_stored = time_temp - time_stored;
            float time_delta = ((float) time_stored)/CLOCKS_PER_SEC;
            printf("Processed %ld reads (%.2f seconds). Estimate time for %.1e reads (%.2f seconds). \n",
                    print_count, time_delta, fix_counter, (time_delta*time_scaler));
            fflush(stdout);

            time_stored = time_temp;

        }
    }


    void CreatePileupV1(boost::program_options::variables_map &vm, GenomeData &genome_data,
                        RefVector &references, SamHeader &header,
                        LocalBamToolsUtils::Fasta &reference_genome,
                        LocalBamToolsUtils::PileupEngine &pileup) {
        BamAlignment ali;
        SampleMap samples;

        uint16_t sindex = 0;
        for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
            if(it->HasSample()){
                auto s  = samples.find(it->Sample);
                if( s == samples.end()){ // not in there yet
                    samples[it->Sample] = sindex;
                    sindex += 1;
                }
            }
        }

        VariantVisitor *v = new VariantVisitor(
                references,
                header,
                reference_genome,
                genome_data,
//            &result_stream,
                samples,
//            params,
                ali,
                vm["qual"].as<int>(),
                vm["mapping-qual"].as<int>(),
                vm["prob"].as<double>()
        );

        pileup.AddVisitor(v);

    }


    void CreatePileupV2(boost::program_options::variables_map &vm, GenomeData &genome_data,
                        RefVector &references, SamHeader &header,
                        LocalBamToolsUtils::Fasta &reference_genome,
                        LocalBamToolsUtils::PileupEngine &pileup) {

        std::string binary_outfile = vm["output_binary_file"].as<std::string>();

        VariantVisitorTwo *v = new VariantVisitorTwo(references, header, reference_genome,
                genome_data, binary_outfile,
                vm["qual"].as<int>(), vm["mapping-qual"].as<int>(), vm["prob"].as<double>());

        pileup.AddVisitor(v);

    }

    void WriteGenomeDataToBinary(std::string file_name, GenomeData &base_counts) {

        GenomeDataStream gd_stream = GenomeDataStream(file_name, true);
//    uint64_t total_base_count = base_counts.size();
        uint64_t sequence_count = base_counts[0].all_reads.size();

        gd_stream.WriteHeader(sequence_count);

        for (auto &baseCount : base_counts) {
            gd_stream.WriteModelInput(baseCount);
        }
        gd_stream.close();
    }

    void ReadGenomeDataFromBinary(std::string file_name, GenomeData &genome_data) {

        GenomeDataStream gd_stream_read = GenomeDataStream( file_name, false);

        uint64_t total_base_count2 = 0;
        uint64_t sequence_count2 = 0;
        gd_stream_read.ReadHeader(total_base_count2, sequence_count2);

        gd_stream_read.ReadGenomeData(genome_data);

        std::cout << "========= Sequence conut: " << sequence_count2  << "\tSite count:" << genome_data.size() << std::endl;
        gd_stream_read.close();

    }



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
//            SequenceProbV1::printReadData(base_counts[i].all_reads[j]) ;

            }
            double prop = (double) sum / total;
            if (prop > 0.1) {
                std::cout << "======= Site: " << i << " R: " << base_counts[i].reference;// << endl;

                std::cout << "\t==" << sum << " " << total << " " << prop << std::endl;

                for (size_t j = 0; j < base_counts[i].all_reads.size(); ++j) {
//                SequenceProbV1::printReadData(base_counts[i].all_reads[j]);
                }
            }
//        cout << "================================="<< endl;
        }
//    cout << "================Done: SequenceProbV1. Total: " << site_count << endl;


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


//    void ReadGenomeDataFromBinaryInit(std::string file_name) {
//
//        GenomeDataStream gd_stream_read = GenomeDataStream( file_name, false);
//
//        uint64_t total_base_count2 = 0;
//        uint64_t sequence_count2 = 0;
//        gd_stream_read.ReadHeader(total_base_count2, sequence_count2);
//        std::cout << "Sequence count: " << sequence_count2 << std::endl;
//
//        gd_stream_read.ReadGenomeData();
//
////        std::cout << "========= conut: " << genome_data.size() << std::endl;
//        gd_stream_read.close();
//
//    }
//
//    void ReadOneRead(){
//        ModelInput m (read_data_count)
//    }

}

