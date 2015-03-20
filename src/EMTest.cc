
#include <iostream>
#include <map>
#include <vector>
#include <boost/program_options.hpp>


#include "api/BamReader.h"
#include "evolution_models/F81.h"
#include "boost_input_utils.h"
#include "pileup_utils.h"
#include "algorithm/em_algorithm_mutation.h"
#include "algorithm/em_algorithm_mutation_v1.h"

using namespace std;
using namespace BamTools;


namespace po = boost::program_options;
int RunBasicProbCalc(GenomeData base_counts, ModelParams params);

void testCalWeighting(MutationProb &mutation_prob, std::vector<SequenceProb> sp);

void AddSimulatedData(ModelParams &params, std::vector<SequenceProb> &sp, int descendant_count, size_t fake_sample_count, double fake_prop);

void SummariseReadsData(GenomeData base_counts);


void RunEmWithRealData(GenomeData &base_counts, ModelParams params) {

    MutationProb mutation_prob = MutationProb(params);
    F81 evo_model0(mutation_prob);

    cout << "init: site_count: " << base_counts.size() << endl;

    std::vector<SequenceProb> sp;
    for (size_t i = 0; i < base_counts.size(); ++i){
    	SequenceProb ss = SequenceProb(base_counts[i], params);
        sp.push_back(ss);
    }

    int descendant_count = sp[0].GetDescendantCount();
    size_t fake_sample_count = 68679;//68679
    double fake_prop = 0.2;//0.2;
    clock_t t1 = clock();
    AddSimulatedData(params, sp, descendant_count, fake_sample_count, fake_prop);
    cout << ((clock() - t1)/ CLOCKS_PER_SEC) << endl;

    cout << "Done preprocess. Final site count: " << sp.size() << endl;

    cout << "======================== Setup EmData:" << endl;



    MutationModel mutation_model = MutationModel(evo_model0);
    mutation_model.AddSequenceProb(sp);
    std::vector<std::unique_ptr<EmModel>> em_model2;
    em_model2.emplace_back(new EmModelMutation(mutation_model));
    em_model2.emplace_back(new EmModelMutation(mutation_model));

    cout << "\n========================\nStart em_algorithm:" << endl;
    clock_t t_start, t_end;

    t_start = clock();
    EmAlgorithmMutation em_alg0 (em_model2);
    em_alg0.Run();
    em_alg0.PrintSummary();
    t_end = clock();
    cout << "Time new: " << (t_end - t_start)/ CLOCKS_PER_SEC << "\t" << (t_end - t_start)  << endl << endl;
}


void AddSimulatedData(ModelParams &params, std::vector<SequenceProb> &sp, int descendant_count, size_t fake_sample_count, double fake_prop) {
    cout << "========= Add simulated data:" << fake_sample_count << endl;

    random_device rd;
    mt19937 e2(rd());
    uniform_int_distribution<uint16_t> uniform_dist(0, 5);
    uniform_int_distribution<uint16_t> uniform3(0, 3);

    size_t fake_diff_count = fake_sample_count * fake_prop;
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



int main(int argc, char** argv){
//using namespace BoostUtils;
    boost::program_options::variables_map variable_map;
    BoostUtils::ParseCommandLinkeInput(argc, argv, variable_map);
//
    GenomeData genome_data;
////    TestVariantVisitorV1(ref_file, vm);
//
////    TestVariantVisitorV1(ref_file, vm, base_counts);
//    PileupUtils::CreatePileupAlignment(variable_map, genome_data, 1);
//    cout << "G0: " << global_count[0] << "\tG1: " << global_count[1] << "\tG2: " << global_count[2] << endl;
//
////    TestVariantVisitorTwo(variable_map, genome_data);
//    PileupUtils::CreatePileupAlignment(variable_map, genome_data, 2);
//    cout << "G0: " << global_count[0] << "\tG1: " << global_count[1] << "\tG2: " << global_count[2] << endl;
//
//
    clock_t start0 = clock();
    std::string file_name = "zz_test_GenomeData_binary_subset";
//    PileupUtils::WriteGenomeDataToBinary(file_name, genome_data);
    PileupUtils::ReadGenomeDataFromBinary(file_name, genome_data);

    printf("Total ReadWrite Time: %f \n", ((double) (clock()-start0)/ CLOCKS_PER_SEC) );
//    SummariseReadsData(genome_data);


    ModelParams params = {
            variable_map["theta"].as<double>(),
            variable_map["nfreqs"].as<vector< double> >(),
            variable_map["mu"].as<double>(),
            variable_map["seq-error"].as<double>(),
            variable_map["phi-haploid"].as<double>(),
            variable_map["phi-diploid"].as<double>(),
    };
//    RunBasicProbCalc(genome_data, params);

    clock_t start = clock();

    RunEmWithRealData(genome_data, params);

    printf("Total EM Time: %f \n", ((double) (clock()-start)/ CLOCKS_PER_SEC) );


}

////////////////////////////////////////////

void SummariseReadsData(GenomeData base_counts){
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


//    EmDataMutationV1 ed2 = EmDataMutationV1();
//    EmData *ed =  &ed2;
//    EmData &ed0 =  ed2;
//
//    EmSummaryStat es;
//    es.stat = 100;
//
//    EmSummaryStatMutationV1 es2;
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

//        std::unique_ptr<EmData> e ( new EmDataMutationV1(sp[s], model) );
//        unique_ptr<EmData> e = std::unique_ptr<EmData>  ( );
//        em_site_data.emplace_back(  new EmDataMutationV1(sp[s], model)  );
//        EmDataMutationV1 em_site = EmDataMutationV1(sp[s], model);
////        EmData *p = &em_site;
//        em_site_prob.push_back(em_site);

    }

    em_site_data.emplace_back(  new EmDataMutationV1(sp[0], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[1], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[2], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[3], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[4], model)  );

    em_site_data.emplace_back(  new EmDataMutationV1(sp[0], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[1], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[2], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[3], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[4], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[0], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[1], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[2], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[3], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[4], model)  );

//    em_site_data.emplace_back(  new EmDataMutationV1(sp[95], model)  );
//    em_site_data.emplace_back(  new EmDataMutationV1(sp[96], model)  );
//    em_site_data.emplace_back(  new EmDataMutationV1(sp[97], model)  );
//    em_site_data.emplace_back(  new EmDataMutationV1(sp[98], model)  );
//    em_site_data.emplace_back(  new EmDataMutationV1(sp[99], model)  );
//    em_site_data.emplace_back(  new EmDataMutationV1(sp[95], model)  );
//    em_site_data.emplace_back(  new EmDataMutationV1(sp[96], model)  );
//    em_site_data.emplace_back(  new EmDataMutationV1(sp[97], model)  );
//    em_site_data.emplace_back(  new EmDataMutationV1(sp[98], model)  );
//    em_site_data.emplace_back(  new EmDataMutationV1(sp[99], model)  );

    em_site_data.emplace_back(  new EmDataMutationV1(sp[95], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[96], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[97], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[98], model)  );
    em_site_data.emplace_back(  new EmDataMutationV1(sp[99], model)  );

    /*

    //Doesn't work, need new
    EmDataMutationV1 em_site = EmDataMutationV1(sp[0], model);
    em_site_prob[0] =(&em_site);
    em_site = EmDataMutationV1(sp[1], model);
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
    EmDataMutationV1 *eeem = static_cast<EmDataMutationV1 *>(e.get());
    eeem->site.CalculateAncestorToDescendant(a,b,c);


    unique_ptr<EmData> &e2 = em_site_data[1];

    EmDataMutationV1 *eeem2 = static_cast<EmDataMutationV1 *>(e2.get());
    eeem2->site.CalculateAncestorToDescendant(a,b,c);



//    exit(100);
    F81 evo_model0(mutation_prob);
//    F81 evo_model1(mutation_prob);
    F81 evo_model1 = evo_model0;
    EmModelMutationV1 em_model0 (evo_model0);
    EmModelMutationV1 em_model1 (em_model0);
//    EmModelMutationV1 em_model1 (evo_model1);

    em_model0.GetParameterInfo();
    em_model1.GetParameterInfo();

//    em_model1.UpdateParameter(0.1);
//    em_model0.UpdateParameter(0.01);
    std::vector<std::unique_ptr<EmModel>> em_model;
    em_model.emplace_back(new EmModelMutationV1(evo_model0));
    em_model.emplace_back(new EmModelMutationV1(em_model0));

//    std::vector<EmModel*> em_model;
//    EmModelMutationV1 em_model1 (em_model0);
//    em_model.push_back(&em_model0);
//    em_model.push_back(&em_model1);

//    EM em (2, site_prob, model);
//    em.RunOld();

//    EmModel *em_base_model = &em_model;
//    EM em (2, site_prob, model, em_model);
//    EM em2 (2, site_prob, model, *em_base_model);

//    EmModel *em_base_model = &em_model;
//    EM em (2, site_prob, model, em_site_prob, em_model);
//    EmModelMutationV1 em_model2 = em_model;

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


//    EmAlgorithmMutationV1 em_alg (2, site_prob, model, em_site_data, em_model0);
    EmAlgorithmMutationV1 em_alg2 (2, em_site_data, em_model0);
    em_alg2.Run();

    em_alg2.PrintSummary();
//    EmAlgorithmMutationV1 em_alg3 (em_site_data, em_model);
//    em_alg3.Run2();

//
////    em_alg.Run();



exit(4);
//    std::vector<std::unique_ptr<EmData>> em_data_mutation;
//    for (size_t s = 0; s < 5; ++s) {
//        em_data_mutation.emplace_back(  new EmDataBinomial(10, s+1)  );
//    }
//    em_data_mutation.emplace_back(  new EmDataBinomial(10, 9)  );
//    em_data_mutation.emplace_back(  new EmDataBinomial(10, 9)  );
//    em_data_mutation.emplace_back(  new EmDataBinomial(10, 10)  );
//    EmModelBinomial em_model_binomial (10, 0.5);
//    EmAlgorithmBinomial em_bin (2, em_data_mutation, em_model_binomial);
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





void TestingEMWithRandomThings(GenomeData base_counts, ModelParams params) {

    MutationProb mutation_prob = MutationProb(params);
    F81 evo_model0(mutation_prob);


    size_t site_count = base_counts.size();
    std::vector<SequenceProb> sp;
//    site_count = 5000;
    cout << "init: site_count: " << site_count << endl;
    for (size_t i = 0; i < site_count; ++i) {
        SequenceProb ss = SequenceProb(base_counts[i], params);
        sp.push_back(ss);
//        sp.push_back(ss);
//        sp.push_back(ss);
    }
//    int descendant_count = sp[0].GetDescendantCount();
//    size_t fake_sample_count = 10000;
//    double fake_prop = 0.25;
//    AddSimulatedData(params, sp, descendant_count, fake_sample_count, fake_prop);

//    TimeTrialWithCache(sp, params);//TODO: Time trial, cache is slower!?!
//    TimeTrialWithCache2(sp, params);//

    cout << "Done preprocess. Final site count: " << sp.size() << endl;


//    MutationModel model = MutationModel(evo_model0);
//    PreprocessSequenceProb2(sp, params, evo_model0);
//    model.AddSequenceProb(sp);
//    double prob, diff, total1, total2;
//    model.cache_read_data_to_all = all_cache_read_to_all_2;
//    for (size_t i = 0; i < sp.size(); ++i) {
//        for (int a = 0; a < 10; ++a) {
//            model.CacheLoopDesAll(i, a, prob, diff);
//            total1 += prob;
//            total2 += diff;
//        }
//        model.CalculateAncestorToDescendant(i, prob, diff);
//        cout << prob << "\t" << diff<< endl;
//    }
//    cout << total1 << "\t" << total2 << endl;
//    exit(11);
//
//    MutationModel model2 = MutationModel(evo_model0);
//    model2.AddSequenceProb(sp);
//    total1 = 0; total2 = 0;
//    for (size_t i = 0; i < sp.size(); ++i) {
//        for (int a = 0; a < 10; ++a) {
//            model2.CacheLoopDesAll(i, a, prob, diff);
//            total1 += prob;
//            total2 += diff;
//        }
//    }
//    cout << total1 << "\t" << total2 << endl;

    int descendant_count = sp[0].GetDescendantCount();
    size_t fake_sample_count = 68679;//68679
    double fake_prop = 0.2;//0.2;
//    AddSimulatedData(params, sp, descendant_count, fake_sample_count, fake_prop);

    cout << "======================== Setup EmData:" << endl;

//    std::vector<SiteProb> site_prob;
    std::vector<EmData*> em_site_prob;
    std::vector<std::unique_ptr<EmData>> em_site_data;
//    for (size_t s = 0; s < site_count; ++s) {
    for (auto seq_prob: sp) {
        em_site_data.emplace_back(  new EmDataMutationV1(seq_prob, evo_model0)  );
//        SiteProb site  (seq_prob, evo_model0 );
//        site_prob.push_back(site);
    }

//    F81 evo_model0(mutation_prob);
//    F81 evo_model1(mutation_prob);
//    F81 evo_model1 = evo_model0;
    EmModelMutationV1 em_model0 (evo_model0);
//    EmModelMutationV1 em_model1 (em_model0);
//    EmModelMutationV1 em_model1 (evo_model1);

//    em_model0.GetParameterInfo();
//    em_model1.GetParameterInfo();

//    em_model1.UpdateParameter(0.1);
//    em_model0.UpdateParameter(0.01);
//    std::vector<std::unique_ptr<EmModel>> em_model;
//    em_model.emplace_back(new EmModelMutationV1(evo_model0));
//    em_model.emplace_back(new EmModelMutationV1(em_model0));


    MutationModel mutation_model = MutationModel(evo_model0);
    mutation_model.AddSequenceProb(sp);
    std::vector<std::unique_ptr<EmModel>> em_model2;
    em_model2.emplace_back(new EmModelMutation(mutation_model));
    em_model2.emplace_back(new EmModelMutation(mutation_model));
//    MutationModel

    cout << "\n========================\nStart em_algorithm:" << endl;
    clock_t t_start, t_end;

    t_start = clock();
    EmAlgorithmMutation em_alg0 (em_model2);
    em_alg0.Run();
    em_alg0.PrintSummary();
    t_end = clock();
    cout << "Time new: " << (t_end - t_start)/ CLOCKS_PER_SEC << endl;

    t_start = clock();
//    EmAlgorithmMutationV1 em_alg (2, site_prob, model, em_site_data, em_model0);
    EmAlgorithmMutationV1 em_alg2 (2, em_site_data, em_model0);
    em_alg2.Run();
    em_alg2.PrintSummary();
    t_end = clock();
    cout << "Time old: " << (t_end - t_start)/CLOCKS_PER_SEC << endl;

//    EmAlgorithmMutationV1 em_alg3 (em_site_data, em_model);
//    em_alg3.Run2();

//

}


