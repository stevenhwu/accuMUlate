/*
 * pileup_utils.cc
 *
 *  Created on: 12/13/14
 *      Author: Steven Wu
 */


#include "pileup_utils.h"

namespace PileupUtils{

    void CreatePileupAlignment(boost::program_options::variables_map &vm, GenomeData &genome_data,
            int variant_visitor_index) {

        BamReader experiment;
        RefVector references;
        SamHeader header;
        Fasta reference_genome;

        BoostUtils::ExtractInputVariables(vm, genome_data, experiment, references, header, reference_genome);
        PileupEngine pileup;

        if(variant_visitor_index == 1){
            CreatePileupV1(vm, genome_data, references, header, reference_genome, pileup);
        }
        else if(variant_visitor_index == 2) {
            CreatePileupV2(vm, genome_data, references, header, reference_genome, pileup);
        }

//  TODO: Only allocate interval-sized memory vector
//  if intervals are set
        BamAlignment ali;
        std::cerr.setstate(std::ios_base::failbit) ; //Supress bamtool camplain, "Pileup::Run() : Data not sorted correctly!"
        if (vm.count("intervals")){
            BedFile bed (vm["intervals"].as<std::string>());
            BedInterval region;
            while(bed.get_interval(region) == 0){
                int ref_id = experiment.GetReferenceID(region.chr);
                experiment.SetRegion(ref_id, region.start, ref_id, region.end);
                while( experiment.GetNextAlignment(ali) ){
                    pileup.AddAlignment(ali);
                }
            }
        }
        else{
            clock_t t;
            clock_t time_temp;
            uint64_t ali_counter = 0;

//        streambuf *old = cout.rdbuf(); // <-- save
//        stringstream ss;
//
//        cout.rdbuf (ss.rdbuf());       // <-- redirect
//        foobar();                      // <-- call
//        cout.rdbuf (old);

            clock_t start = clock();
            t=start;
            global_count[0] = 0;global_count[1] = 0;global_count[2] = 0;

            while( experiment.GetNextAlignment(ali)){           // Fast, 0.2s
                pileup.AddAlignment(ali);                     // AddAlignment ~2s + visitor ~3s
                ali_counter += 1;
                if (ali_counter % 100000 == 0){  //3800 => 67
                time_temp = clock() ;
                t = time_temp - t;
                    std::cout << "Processed 1 million reads ("
                        << ((float)t)/CLOCKS_PER_SEC
                        << " seconds)" << std::endl;
                t = time_temp;
//                    break;
                }

            }
            pileup.Flush();

            printf("Total time V%d: %f Count:%lu G_count:%d Base_count:%lu\n", variant_visitor_index,
                    ((double) (clock()-start)/ CLOCKS_PER_SEC),
                    ali_counter, global_count[0], genome_data.size());
            std::cout << "====DEBUG: G0: " << global_count[0] << "\tG1: " << global_count[1] << "\tG2: " << global_count[2] << std::endl;
            std::cout.flush();

/*
        Let's work with 50k count for now, 4.5~5s   6788 base_count, 2s overhead
        //Init test: ~1 to 1.5s per 10k.
        // ~25s for 30.7k base_count, 194927 ali_count
        // 10s for 194927 count without visitor
        update: 50k ~ 3.3~4s
        Full 194927 ali_count: 15-18s
        Full EM: ~130s

*/
        }
        std::cerr.clear() ;
        experiment.Close();
        reference_genome.Close();
    }


    void CreatePileupV1(boost::program_options::variables_map &vm, GenomeData &genome_data, RefVector &references, SamHeader &header, Fasta &reference_genome, PileupEngine &pileup) {
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


    void CreatePileupV2(boost::program_options::variables_map &vm, GenomeData &genome_data, RefVector &references, SamHeader &header, Fasta &reference_genome, PileupEngine &pileup) {
        VariantVisitorTwo *v = new VariantVisitorTwo(
                references,
                header,
                reference_genome,
                genome_data,
//            &result_stream,
//            samples,
//            params,
//            ali,
                vm["qual"].as<int>(),
                vm["mapping-qual"].as<int>(),
                vm["prob"].as<double>()
        );

        pileup.AddVisitor(v);
    }


    void WriteGenomeDataToBinary(std::string file_name, GenomeData &base_counts) {

        GenomeDataStream gd_stream = GenomeDataStream(file_name, true);
//    uint64_t total_base_count = base_counts.size();
        uint64_t sequence_count = base_counts[0].all_reads.size();

        gd_stream.WriteHeader(sequence_count);

        for (auto baseCount : base_counts) {
            gd_stream.WriteModelInput(baseCount);
        }
        gd_stream.close();
    }

    void ReadGenomeDataFromBinary(std::string file_name, GenomeData &genome_data) {

        GenomeDataStream gd_stream_read = GenomeDataStream( file_name, false);

        uint64_t total_base_count2 = 0;
        uint64_t sequence_count2 = 0;
        gd_stream_read.ReadHeader(total_base_count2, sequence_count2);
        std::cout << total_base_count2 << "\t" << sequence_count2 << std::endl;

        gd_stream_read.ReadGenomeData(genome_data);


//    int count = 0;
//    for (auto baseCount : base_counts) {
////        cout << baseCount.reference << " : ";
//        for (auto item : baseCount.all_reads) {
////            cout << item.key << " : ";
////            SequenceProb::printReadData(item);
//        }
////        cout << "\n";
//        count++;
////        if(count == 2){
//////            break;
////        }
//    }

        std::cout << "========= conut: " << genome_data.size() << std::endl;
        gd_stream_read.close();

    }




}