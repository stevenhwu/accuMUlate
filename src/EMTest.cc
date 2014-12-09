
#include <iostream>
#include <map>
#include <vector>

#include <boost/program_options.hpp>
#include <time.h>
#include <algorithm/em_model_binomial.h>
#include <algorithm/em_data_binomial.h>
#include <algorithm/em_algorithm_binomial.h>


#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"
#include "parsers.h"
#include "sequence_prob.h"
#include "VariantVisitor.h"
#include "evolution_models/JC69.h"
#include "evolution_models/F81.h"
#include "site_prob.h"
#include "algorithm/em_algorithm_mutation.h"
#include "algorithm/em_data_mutation.h"

using namespace std;
using namespace BamTools;


namespace po = boost::program_options;
int RunBasicProbCalc(GenomeData base_counts, ModelParams params);

void testCalLikelihood(MutationProb muProb, std::vector<SequenceProb> sp);

void testCalWeighting(MutationProb &mutation_prob, std::vector<SequenceProb> sp);



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
//    F81  model(muArray[0], model_params.nuc_freq);
    site_count = 100;
    for (size_t i = 0; i < site_count; ++i) {
        sp.push_back(SequenceProb(base_counts[500 + i], params));
    }
//    sp[0] = SequenceProb(base_counts[598], params);
//    sp[1] = SequenceProb(base_counts[4195], params);
//		sp[2] = SequenceProb(params, base_counts[5391], muProb);
//		sp[3] = SequenceProb(params, base_counts[6589], muProb);

    ModelInput base_custom =  base_counts[598];

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1,10);
    int dice_roll = distribution(generator);

    std::random_device rd;
//    std::default_random_engine e1(rd());
//    std::uniform_int_distribution<int> uniform_dist(1, 6);
//    int mean = uniform_dist(e1);

    // Generate a normal distribution around that mean
    std::mt19937 e2(rd());
    std::uniform_int_distribution<int> uniform_dist(5, 10);
    for (size_t s = 0; s < site_count; ++s) {
        ModelInput base_custom = base_counts[s];
        base_custom.reference = 0;
        for (int i = 0; i < 7; ++i) {
            base_custom.all_reads[i].key=0;
            for (int j = 0; j < 4; ++j) {
                base_custom.all_reads[i].reads[j] = (uint16_t) uniform_dist(e2);
            }
            base_custom.all_reads[i].reads[0] = (uint16_t) 100 + uniform_dist(e2);
        }
//        int i=0;
        if(s > 70){ //diff
            base_custom.all_reads[0].reads[0] = (uint16_t) uniform_dist(e2);
            base_custom.all_reads[0].reads[3] = (uint16_t) 100 + uniform_dist(e2);
            base_custom.all_reads[1].reads[2] = (uint16_t) 100 + uniform_dist(e2);
//            base_custom.all_reads[2].reads[3] = (uint16_t) 100 + uniform_dist(e2);
//            base_custom.all_reads[3].reads[0] = (uint16_t) 100 + uniform_dist(e2);
//            base_custom.all_reads[4].reads[1] = (uint16_t) 100 + uniform_dist(e2);
//            base_custom.all_reads[5].reads[2] = (uint16_t) 100 + uniform_dist(e2);
//            base_custom.all_reads[6].reads[3] = (uint16_t) 100 + uniform_dist(e2);
//            base_custom.all_reads[0].reads[2] = (uint16_t) 100 + uniform_dist(e2);
        }
        sp[s] = SequenceProb(base_custom, params);
    }

//    for (int s = 75; s < site_count; ++s) {
//        ModelInput base_custom = base_counts[s];
//        base_custom.reference = 0;
//        for (int i = 0; i < 7; ++i) {
//            base_custom.all_reads[i].key=0;
//            for (int j = 0; j < 4; ++j) {
//                base_custom.all_reads[i].reads[j] = (uint16_t) uniform_dist(e2);
//            }
//            base_custom.all_reads[i].reads[2] = (uint16_t) 100+uniform_dist(e2);
//        }
//
//
//
//        sp[s] = SequenceProb(base_custom, params);
//    }


    cout << "Start\n\n";


//    EmDataMutation ed2 = EmDataMutation();
//    EmData *ed =  &ed2;
//    EmData &ed0 =  ed2;
//
//    EmSummaryStat es;
//    es.stat = 100;
//
//    EmSummaryStatMutation es2;
//    es2.stat = 200;
//    es2.stat_diff = 210;
//
//
//    ed->UpdateSummaryStat(0, es);
//    ed2.UpdateSummaryStat(0, es);
//    ed0.UpdateSummaryStat(0, es);
//
//    ed->UpdateSummaryStat(0, es2);
//    ed2.UpdateSummaryStat(0, es2);
//    ed0.UpdateSummaryStat(0, es2);
//
//    es2.print();


    testCalWeighting(muProb, sp);

    return 0;
}

void testCalWeighting(MutationProb &mutation_prob, std::vector<SequenceProb> sp) {
    const size_t cat = 2;


    double muArray[cat];
    muArray[0] = 1e-10;
    muArray[1] = 1;

//    double muNewArray[cat];
//    muNewArray[0] = muArray[0];
//    muNewArray[1] = muArray[1];


    double proportion[cat] {0.5,0.5};
    double all_stats_same[cat];
    double all_stats_diff[cat];

//    double weight[2];


//    JC69 model(muArray[0]);
//    MutationProb muProbs[cat] = {
//            MutationProb(model_params),
//            MutationProb(model_params)};
    //F81 model2(muArray[0], model_params.nuc_freq);

//    F81 model2(0.1);
//    JC69 m69(0.1);
//
//    JC69 m692(mutation_prob);

    F81 model(mutation_prob);

//    exit(-10);
//    F81 model(mutation_prob);
//    JC69 m69(mutation_prob);
//    MutationMatrix conditional_prob = model.GetTransitionMatirxAToD();
    const int c_site_count = 100;
    size_t site_count = c_site_count;
    std::vector<SiteProb> site_prob;
    std::vector<EmData*> em_site_prob;

    std::vector<std::unique_ptr<EmData>> em_site_data;
//    stuff.reserve(site_count);
//    for( int i = 0; i < 10; ++i )


    for (size_t s = 0; s < site_count; ++s) {
//        SiteProb site  (sp[s],mutation_prob, model );
        SiteProb site  (sp[s], model );
        site_prob.push_back(site);

//        std::unique_ptr<EmData> e ( new EmDataMutation(sp[s], model) );
//        unique_ptr<EmData> e = std::unique_ptr<EmData>  ( );
//        em_site_data.emplace_back(  new EmDataMutation(sp[s], model)  );
//        EmDataMutation em_site = EmDataMutation(sp[s], model);
////        EmData *p = &em_site;
//        em_site_prob.push_back(em_site);

    }

    em_site_data.emplace_back(  new EmDataMutation(sp[0], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[1], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[2], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[3], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[4], model)  );

    em_site_data.emplace_back(  new EmDataMutation(sp[0], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[1], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[2], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[3], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[4], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[0], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[1], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[2], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[3], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[4], model)  );

//    em_site_data.emplace_back(  new EmDataMutation(sp[95], model)  );
//    em_site_data.emplace_back(  new EmDataMutation(sp[96], model)  );
//    em_site_data.emplace_back(  new EmDataMutation(sp[97], model)  );
//    em_site_data.emplace_back(  new EmDataMutation(sp[98], model)  );
//    em_site_data.emplace_back(  new EmDataMutation(sp[99], model)  );
//    em_site_data.emplace_back(  new EmDataMutation(sp[95], model)  );
//    em_site_data.emplace_back(  new EmDataMutation(sp[96], model)  );
//    em_site_data.emplace_back(  new EmDataMutation(sp[97], model)  );
//    em_site_data.emplace_back(  new EmDataMutation(sp[98], model)  );
//    em_site_data.emplace_back(  new EmDataMutation(sp[99], model)  );

    em_site_data.emplace_back(  new EmDataMutation(sp[95], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[96], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[97], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[98], model)  );
    em_site_data.emplace_back(  new EmDataMutation(sp[99], model)  );

    /*

    //Doesn't work, need new
    EmDataMutation em_site = EmDataMutation(sp[0], model);
    em_site_prob[0] =(&em_site);
    em_site = EmDataMutation(sp[1], model);
    em_site_prob[1] = (&em_site);
*/
    size_t em_count = 50;
    size_t rate_count = 2; //2;
    double all_prob[cat][c_site_count];


    cout << "\n========================\nStart em_algorithm:" << endl;
    double a,b,c;
    cout << site_prob.size() << "\t" << em_site_prob.size() << endl;

    site_prob[0].CalculateAncestorToDescendant(a,b,c);
    site_prob[1].CalculateAncestorToDescendant(a,b,c);
    cout <<"\n\n";
    unique_ptr<EmData> &e = em_site_data[0];
//    EmData *ed= e.get();
    EmDataMutation *eeem = static_cast<EmDataMutation*>(e.get());
    eeem->site.CalculateAncestorToDescendant(a,b,c);


    unique_ptr<EmData> &e2 = em_site_data[1];

    EmDataMutation *eeem2 = static_cast<EmDataMutation*>(e2.get());
    eeem2->site.CalculateAncestorToDescendant(a,b,c);



//    exit(100);
    F81 evo_model0(mutation_prob);
//    F81 evo_model1(mutation_prob);
    F81 evo_model1 = evo_model0;
    EmModelMutation em_model0 (evo_model0);
    EmModelMutation em_model1 (em_model0);
//    EmModelMutation em_model1 (evo_model1);

    em_model0.GetParameterInfo();
    em_model1.GetParameterInfo();

//    em_model1.UpdateParameter(0.1);
//    em_model0.UpdateParameter(0.01);
    std::vector<std::unique_ptr<EmModel>> em_model;
    em_model.emplace_back(new EmModelMutation(evo_model0));
    em_model.emplace_back(new EmModelMutation(em_model0));

//    std::vector<EmModel*> em_model;
//    EmModelMutation em_model1 (em_model0);
//    em_model.push_back(&em_model0);
//    em_model.push_back(&em_model1);

//    EM em (2, site_prob, model);
//    em.RunOld();

//    EmModel *em_base_model = &em_model;
//    EM em (2, site_prob, model, em_model);
//    EM em2 (2, site_prob, model, *em_base_model);

//    EmModel *em_base_model = &em_model;
//    EM em (2, site_prob, model, em_site_prob, em_model);
//    EmModelMutation em_model2 = em_model;

    em_model[0]->GetParameterInfo();
    em_model[1]->GetParameterInfo();
    printf("H1\n");

    em_model[1]->UpdateParameter(0.1);

    em_model[0]->GetParameterInfo();
    em_model[1]->GetParameterInfo();
    printf("H2\n");

    em_model[0]->UpdateParameter(1e-2);

    em_model[0]->GetParameterInfo();
    em_model[1]->GetParameterInfo();
    printf("H3\n");


//    EmAlgorithmMutation em_alg (2, site_prob, model, em_site_data, em_model0);
    EmAlgorithmMutation em_alg2 (2, em_site_data, em_model0);
    em_alg2.Run();

//    EmAlgorithmMutation em_alg3 (em_site_data, em_model);
//    em_alg3.Run2();

//
////    em_alg.Run();



exit(4);
//    std::vector<std::unique_ptr<EmData>> em_data_binomial;
//    for (size_t s = 0; s < 5; ++s) {
//        em_data_binomial.emplace_back(  new EmDataBinomial(10, s+1)  );
//    }
//    em_data_binomial.emplace_back(  new EmDataBinomial(10, 9)  );
//    em_data_binomial.emplace_back(  new EmDataBinomial(10, 9)  );
//    em_data_binomial.emplace_back(  new EmDataBinomial(10, 10)  );
//    EmModelBinomial em_model_binomial (10, 0.5);
//    EmAlgorithmBinomial em_bin (2, em_data_binomial, em_model_binomial);
//
//
//    em_bin.Run();
//
//    vector<double> p1 = em_bin.GetParameters();
//    for (auto item : p1) {
//        printf("%f\t", item);
//    }
//    printf("\n");

//exit(35);


    for (size_t i = 0; i < em_count; ++i) {

        for (size_t r = 0; r < rate_count; ++r) {
            model.UpdateMu(muArray[r]);
//            mutation_prob.UpdateMu(muArray[r]);

            all_stats_same[r] = 0;
            all_stats_diff[r] = 0;
//            all_prob[r][s] = 0;
//          Array10D site_stat[site_count];
            for (size_t s = 0; s < site_count; ++s) {

                auto site = site_prob[s];
//                t.UpdateTransitionMatrix(model);
//                t.UpdateMuProb(mutation_prob);
                site.UpdateModel(model);
                double stat_same = 0;
                double stat_diff = 0;
                double sum_prob = 0;
                site.CalculateAncestorToDescendant(sum_prob, stat_same, stat_diff);
                all_stats_same[r] += proportion[r] * stat_same;
                all_stats_diff[r] += proportion[r] * stat_diff;
                all_prob[r][s] = sum_prob;

            }
        }

        double sum = 0;
        for (size_t s = 0; s < c_site_count; ++s) {
            double sum0 = all_prob[0][s]/ (all_prob[0][s]+all_prob[1][s]);
            sum += sum0;
        }
        proportion[0] = sum/site_count;
        proportion[1] = 1-proportion[0];

        for (size_t r = 0; r < rate_count; ++r) {
            double sum_stat = all_stats_diff[r] + all_stats_same[r];
            double new_exp_beta = all_stats_diff[r] / sum_stat;
            double new_mu = mutation_prob.ConvertExpBetaToMu(new_exp_beta);

            double new_one_minus_exp_beta = all_stats_same[r] / sum_stat;
            double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta-1));
//            cout.setprecision(10);

//            proportion[r] = all_prob[r]/(all_prob[0]+all_prob[1]);
            printf("======= NEM_MU_r: %zu \tMu: %.5e %.5f =expBeta=  %.5f %.5f \t =prop= %.5f %.5f "
                            "\t =all_prob=  %.5f \t =stat= %.5f %.5f \n" ,r,
                    new_mu ,new_one_minus_mu,
                    new_exp_beta, new_one_minus_exp_beta,
                    proportion[0], proportion[1],  sum,
                    all_stats_same[0], all_stats_same[1]);

            muArray[r] = new_mu;

        }
        

    }
    

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
                SequenceProb sp(base, params);
//                sp.UpdateSummaryStat();
//                double likelihood = sp.GetLikelihood();
//                cout << likelihood << (1 - likelihood) << endl;

                const int size = 20;
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

