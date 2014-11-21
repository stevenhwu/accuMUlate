
#include <iostream>
#include <map>
#include <vector>

#include <boost/program_options.hpp>
#include <models/JC69.h>

#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"
#include "parsers.h"
#include "SequenceProb.h"
#include "VariantVisitor.h"
#include "models/F81.h"

using namespace std;
using namespace BamTools;


namespace po = boost::program_options;
int RunBasicProbCalc(GenomeData base_counts, ModelParams params);

void testCalLikelihood(MutationProb muProb, std::vector<SequenceProb> sp);

void testCalWeighting(MutationProb muProb, std::vector<SequenceProb> sp, ModelParams model_params);



int main(int argc, char** argv){

    namespace po = boost::program_options;
    string ref_file;
    string config_path;
    po::options_description cmd("Command line options");
    cmd.add_options()
            ("help,h", "Print a help message")
            ("bam,b", po::value<string>()->required(), "Path to BAM file")
            ("bam-index,x", po::value<string>()->default_value(""), "Path to BAM index, (defalult is <bam_path>.bai")
            ("reference,r", po::value<string>(&ref_file)->required(),  "Path to reference genome")
//       ("ancestor,a", po::value<string>(&anc_tag), "Ancestor RG sample ID")
//        ("sample-name,s", po::value<vector <string> >()->required(), "Sample tags")
            ("qual,q", po::value<int>()->default_value(13),
                    "Base quality cuttoff")

            ("mapping-qual,m", po::value<int>()->default_value(13),
                    "Mapping quality cuttoff")

            ("prob,p", po::value<double>()->default_value(0.1),
                    "Mutaton probability cut-off")
            ("out,o", po::value<string>()->default_value("acuMUlate_result.tsv"),
                    "Out file name")
            ("intervals,i", po::value<string>(), "Path to bed file")
            ("config,c", po::value<string>(), "Path to config file")
            ("theta", po::value<double>()->required(), "theta")
            ("nfreqs", po::value<vector<double> >()->multitoken(), "")
            ("mu", po::value<double>()->required(), "")
            ("seq-error", po::value<double>()->required(), "")
            ("phi-haploid",     po::value<double>()->required(), "")
            ("phi-diploid",     po::value<double>()->required(), "");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmd), vm);

    if (vm.count("help")){
        cout << cmd << endl;
        return 0;
    }

    if (vm.count("config")){
        ifstream config_stream (vm["config"].as<string>());
        po::store(po::parse_config_file(config_stream, cmd, false), vm);
    }

    vm.notify();
//    ModelParams params = {
//        vm["theta"].as<double>(),
//        vm["nfreqs"].as<vector< double> >(),
//        vm["mu"].as<double>(),
//        vm["seq-error"].as<double>(),
//        vm["phi-haploid"].as<double>(),
//        vm["phi-diploid"].as<double>(),
//    };
    string bam_path = vm["bam"].as<string>();
    string index_path = vm["bam-index"].as<string>();
    if(index_path == ""){
        index_path = bam_path + ".bai";
    }

//    ofstream result_stream (vm["out"].as<string>());
    //TODO: check sucsess of all these opens/reads:
    BamReader experiment;
    experiment.Open(bam_path);
    experiment.OpenIndex(index_path);
    RefVector references = experiment.GetReferenceData();
    SamHeader header = experiment.GetHeader();
    Fasta reference_genome; // BamTools::Fasef_file);
    reference_genome.Open(ref_file, ref_file+ ".fai");
//    reference_genome.CreateIndex(ref_file + ".fai");
    PileupEngine pileup;
    BamAlignment ali;


//  Assign some memory for the big list
    uint32_t total_len = 0;
    if (vm.count("intervals")){
        BedFile bed (vm["intervals"].as<string>());
        BedInterval region;
        while(bed.get_interval(region) == 0){
            total_len += (region.end - region.start);
        }
    }
    else {

        for_each(references.begin(),references.end(),[&](RefData chrom){
            total_len += chrom.RefLength;
        });
    }
    GenomeData base_counts;
    base_counts.reserve(total_len);

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
            base_counts,
//            &result_stream,
            samples,
//            params,
            ali,
            vm["qual"].as<int>(),
            vm["mapping-qual"].as<int>(),
            vm["prob"].as<double>()
    );
    pileup.AddVisitor(v);
//  TODO: Only allocate interval-sized memory vector
//  if intervals are set
    if (vm.count("intervals")){
        BedFile bed (vm["intervals"].as<string>());
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
        uint64_t ali_counter = 0;
        t = clock();
        BamAlignment ali;

        std::cerr.setstate(std::ios::failbit) ;
//        streambuf *old = cout.rdbuf(); // <-- save
//        stringstream ss;
//
//        cout.rdbuf (ss.rdbuf());       // <-- redirect
//        foobar();                      // <-- call
//        cout.rdbuf (old);

        while( experiment.GetNextAlignment(ali)){
            pileup.AddAlignment(ali);
            ali_counter += 1;
            if (ali_counter % 1000000 == 0){
                t = clock() - t;
                cout << "Processed 1 million reads ("
                        << ((float)t)/CLOCKS_PER_SEC
                        << " seconds)" << endl;
            }
        }
        std::cerr.clear() ;
    }
    pileup.Flush();
    cout << "Base_count: " << base_counts.size() << endl;

    ModelParams params = {
        vm["theta"].as<double>(),
        vm["nfreqs"].as<vector< double> >(),
        vm["mu"].as<double>(),
        vm["seq-error"].as<double>(),
        vm["phi-haploid"].as<double>(),
        vm["phi-diploid"].as<double>(),
    };
    RunBasicProbCalc(base_counts, params);

}

int RunBasicProbCalc(GenomeData base_counts, ModelParams params) {
    cout << "init" << endl;
    size_t site_count = base_counts.size();
    MutationProb muProb = MutationProb(params);
    std::vector<SequenceProb> sp;

//    for (int i = 0; i < site_count; ++i) {
//        sp[i] = SequenceProb(params, base_counts[+i], muProb);
//    }
    for (int i = 0; i < 3; ++i) {
        sp.push_back( SequenceProb(params, base_counts[500+i], muProb) );
    }
    sp[0] = SequenceProb(params, base_counts[598], muProb);
//		sp[1] = SequenceProb(params, base_counts[4195], muProb);
    sp[1] = SequenceProb(params, base_counts[5391], muProb);
//		sp[2] = SequenceProb(params, base_counts[5391], muProb);
//		sp[3] = SequenceProb(params, base_counts[6589], muProb);

    ModelInput base_custom =  base_counts[598];

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1,10);

    std::random_device rd;
//    std::default_random_engine e1(rd());
//    std::uniform_int_distribution<int> uniform_dist(1, 6);
//    int mean = uniform_dist(e1);

    // Generate a normal distribution around that mean
    std::mt19937 e2(rd());
    std::uniform_int_distribution<int> uniform_dist(1, 30);

    int dice_roll = distribution(generator);
    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 4; ++j) {
            base_custom.all_reads[i].reads[j] = (uint16_t) uniform_dist(e2);
        }
//        base_custom.all_reads[i].reads[0] = (uint16_t) i+0;
    }


    sp[2] = SequenceProb(params, base_custom, muProb);

//		double d = sp[0].GetLikelihood();
		cout << "D:\t"  << endl;
		for (auto t : sp) {
			t.UpdateLikelihood();
			double likelihood = t.GetLikelihood();
			cout << likelihood << (1 - likelihood) << endl;
		}

    cout << "Start\n\n";

//    testCalLikelihood(muProb, sp);
    testCalWeighting(muProb, sp, params);

    return 0;
}

void testCalWeighting(MutationProb muProb, std::vector<SequenceProb> sp, ModelParams model_params) {
    const size_t cat = 2;
    size_t site_count = 1000;//

    double muArray[cat];
    muArray[0] = 1e-2;
    muArray[1] = 1e-10;

    double proportion[cat];
    double weight[2];
    int num_descendant = 2;//FIXME change later

//    JC69 model(muArray[0]);
    F81  model(muArray[0], model_params.nuc_freq);
//    MutationMatrix conditional_prob = model.GetTransitionMatirxAToD();
    for (int r = 0; r < 1; ++r) {
        model.UpdateMu(muArray[r]);
        muProb.UpdateMu(muArray[r]);

        Array10D site_stat[site_count];
        for (int s = 0; s < 1; ++s) {
            s=2;
            auto t = sp[s];
            cout << "==============Start Looop:\t Site: " << s << " rate:" << r  << "\n";

            t.UpdateTransitionMatrix(model);
            t.UpdateMuProb(muProb);

//            for (int b = 0; b < 4; ++b) {
                
//                cout << "here1\n";

//                prob_AtoD +=
            t.CalculateAncestorToDescendant();

//                prob_reads_given_descent
//                cout << "here2\n";
//                cout << prob_for_each_base << endl;
//                cout << prob_reads_given_descent << "\t" << prob_reads_given_descent[0]<<
//                "\t" << prob_reads_given_descent[1]<<
//                "\t" << prob_reads_given_descent[2]<<
//                "\t" << prob_reads_given_descent[3]<< endl;

//                likelihood[i] = t.GetLikelihood();
//            }


//            cout << "====================================summary_stat: " << summary_stat  << endl;
        }
//        double sum = likelihood[0] + likelihood[1];
//
//        proportion[0] = likelihood[0]/sum;
//        proportion[1] = likelihood[1]/sum;
//
//        cout << proportion[0] << "\t" << proportion[1] << endl;

//        ModelInput m = t.GetData();
////			m.all_reads
//
//        ModelInput mA, mB;
//        ReadData vA, vB;
//        //Might have to redo this part, with double[4] TODO:
//        mA.all_reads = std::vector<ReadData>(7);
        /*for (int b = 0; b < m.all_reads.size(); ++b) {
            ReadData v = m.all_reads[b];

            cout << v.reads[0] << "\t" << v.reads[1] << "\t"<< v.reads[2] << "\t"<< v.reads[3] << "\t" <<endl;
            for (int bb = 0; bb < 3; ++bb) {
                vA.reads[bb] = v.reads[bb] * proportion[0];
                vB.reads[bb] = v.reads[bb] * proportion[1];
            }
            cout << vA.reads[0] << "\t" << vA.reads[1] << "\t"<< vA.reads[2] << "\t"<< vA.reads[3] << "\t" <<endl;
            cout << "AA" << endl;
            mA.all_reads[b] = vA;
            cout << "BB" << endl;
            mB.all_reads.push_back( vB );
            cout << "??" << endl;
        }*/
    }
//    cout << muArray[0] << "\t" << muArray[2] << endl;
}

void testCalLikelihood(MutationProb muProb, SequenceProb sp[]) {
    int size = 5;
    double muArray[size];
    double likelihood[size];
    muArray[0] = 1;
    for (int i = 1; i < size; ++i) {
        muArray[i] = muArray[i - 1] * 0.1;
    }
    for (int i = 0; i < size; ++i) {
        muProb.UpdateMu(muArray[i]);
        for (int s = 0; s < 3; ++s) {
            auto t = sp[s];
            t.UpdateMuProb(muProb);
            cout << log( t.GetLikelihood() )<< " ";
            likelihood[i] += log( t.GetLikelihood() );
        }
        printf("\tI:%d\t%.1e %e %f \n", i, muArray[i], likelihood[i], likelihood[i]
        );
    }
    double mma = *max_element( likelihood, likelihood+size);
    double mmi = *min_element( likelihood, likelihood+size);
    cout << mma <<"\t"<< mmi << endl;
}

void RunMaProb(ModelParams params, po::variables_map vm, BamReader experiment, PileupEngine pileup, GenomeData base_counts) {
    BamAlignment ali;
    if (vm.count("intervals")){
        BedFile bed (vm["intervals"].as<string>());
        BedInterval region;
        while(bed.get_interval(region) == 0){
            int ref_id = experiment.GetReferenceID(region.chr);
            experiment.SetRegion(ref_id, region.start, ref_id, region.end);
            while( experiment.GetNextAlignment(ali) ){
                pileup.AddAlignment(ali);
            }
        }
    }
    else {
        cout << "here" << endl;


        cout << base_counts.size() << endl;
        cout << "here3" << endl;

//        ModelParams params = {
//                vm["theta"].as<double>(),
//                vm["nfreqs"].as<vector<double> >(),
//                vm["mu"].as<double>(),
//                vm["seq-error"].as<double>(),
//                vm["phi-haploid"].as<double>(),
//                vm["phi-diploid"].as<double>(),
//        };
        int count = 0;
        for (auto base : base_counts) {

//		ModelInput d = {ref_base_idx, bcalls};
//		double prob_one = TetMAProbOneMutation(params,d);
            for (auto t : base.all_reads) {
                cout << t.reads[0] << t.reads[1] << t.reads[2] << t.reads[3] << endl;
            }
            cout << "" << endl;
//    	exit(-1);
//            double prob = TetMAProbability(params, base);
            double prob = 0;
            if (prob >= 0.1) {
                cout.precision(10);
                cout << "=================" << endl;
                cout << prob << "\t" << (1 - prob) << endl;
                cout << base.reference << "\tSIZE:\t" << base.all_reads.size() << endl;
                for (auto t : base.all_reads) {
                    cout << t.reads[0] << t.reads[1] << t.reads[2] << t.reads[3] << endl;
                }
                cout << count << endl;
                MutationProb muProb = MutationProb(params);
                SequenceProb sp(params, base, muProb);
                sp.UpdateLikelihood();
//                double likelihood = sp.GetLikelihood();
//                cout << likelihood << (1 - likelihood) << endl;

                int size = 20;
                double muArray[size];
                muArray[0] = 0.1;
                for (int i = 1; i < size; ++i) {
                    muArray[i] = muArray[i - 1] * 0.1;
                }
                for (int i = 0; i < size; ++i) {
//                    muProb.UpdateMu(muArray[i]);
//                    sp.UpdateMuProb(muProb);
//                    likelihood = sp.GetLikelihood();
//
//                    printf("I:%d %e %e %e\n", i, muArray[i], likelihood, (1 - likelihood));
                }



//			sp.UpdateMu(1e-7);
//			 likelihood = sp.GetLikelihood();
//			cout << likelihood << (1-likelihood) << endl;
//
//			sp.UpdateMu(1e-8);
//			 likelihood = sp.GetLikelihood();
//			cout << likelihood << (1-likelihood) << endl;
//
//			sp.UpdateMu(1e-60);
//			 likelihood = sp.GetLikelihood();
//			cout << likelihood << (1-likelihood) << endl;


                cout << "=================" << endl;
////					std::cout << base.reference << base.all_reads << endl;
////			 *m_ostream << chr << '\t'
////						<< pos << '\t'
////						<< current_base << '\t'
////						<< prob << '\t'
////						<< prob_one << '\t'
////					   << endl;
            }

            count++;
            if (count == 600) { //598 4195 5391 6589
                break;
            }
        }
    }
}

