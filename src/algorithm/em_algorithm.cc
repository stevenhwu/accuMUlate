

#include "em_algorithm.h"
#include "em_data_mutation.h"
#include "em_model_mutation.h"

EM::EM(int num_category0, vector<SiteProb> em_data0, EvolutionModel &em_model0) : num_category(num_category0), em_data(em_data0), em_model(&em_model0){

    cout << "EM Constructor\n";

    if (num_category != 2) {
        cout << "Not yet implemented for more than 2 categories" << endl;
        exit(222);
    }

    site_count = em_data.size();
    rate_count = num_category;
    em_count = 10;

    Init();



}

//EM::EM(vector<EmDataMutation> site_data, EmModelMutation &evo_model){//: em_data(em_data), em_model(&em_model){
//
//    site_count = site_data.size();
//}


EM::~EM() {

}

void EM::Init() {

    InitialiseParameters();
    InitialiseProportion();
    all_probs = Eigen::ArrayXXd::Zero(num_category, site_count);


    all_em_stats = vector<EmSummaryStatMutation>(rate_count);

    all_stats_same = vector<double>(num_category, 0);
    all_stats_diff = vector<double>(num_category, 0);



}


void EM::ExpectationStep() {

    for (size_t r = 0; r < rate_count; ++r) {
        em_model->UpdateMu(parameters[r]);
//            mutation_prob.UpdateMu(parameters[r]);

        all_em_stats[r].Reset();

        for (size_t s = 0; s < site_count; ++s) {
            double sum_prob = 0;

            auto site = em_data[s];

            site.UpdateModel(*em_model);

            double stat_same = 0;
            double stat_diff = 0;
            //relpace with
            EmSummaryStatMutation local;


            site.CalculateAncestorToDescendant(sum_prob, stat_same, stat_diff);


            local.SetStats(stat_same, stat_diff);
            all_em_stats[r].UpdateSumWithProportion(proportion[r], local);


            all_probs(r,s) = sum_prob;

        }

    }


}

void EM::Run() {

    MutationProb mutation_prob = em_model-> GetMutationProb();


//    for (size_t s = 0; s < site_count; ++s) {
//        EmSummaryStatMutation e;
////        e.stat = 1;
////        e.stat_diff = 2;
//        all_em_stats.push_back(e);
//    }

em_count = 50;

    for (size_t i = 0; i < em_count; ++i) {
cout << "EM ite: " << i << endl;

        for (size_t r = 0; r < rate_count; ++r) {
            em_model->UpdateMu(parameters[r]);
//            mutation_prob.UpdateMu(parameters[r]);

            all_stats_same[r] = 0;
            all_stats_diff[r] = 0;
            //replace with
//            auto em_stat = all_em_stats[r];
            all_em_stats[r].Reset();

//            em_stat.print();
//            all_probs[r][s] = 0;
//          Array10D site_stat[site_count];
            for (size_t s = 0; s < site_count; ++s) {
                double sum_prob = 0;

                auto site = em_data[s];

                site.UpdateModel(*em_model);

                double stat_same = 0;
                double stat_diff = 0;
                //relpace with
                EmSummaryStatMutation local;


                site.CalculateAncestorToDescendant(sum_prob, stat_same, stat_diff);


                all_stats_same[r] += proportion[r] * stat_same;
                all_stats_diff[r] += proportion[r] * stat_diff;
                //replace with
                local.SetStats(stat_same, stat_diff);
                all_em_stats[r].UpdateSumWithProportion(proportion[r], local);


                if(all_em_stats[r].stat != all_stats_same[r] || all_em_stats[r].stat_diff!=all_stats_diff[r]){
                    cout << "DIFF!!! inner " << endl;
                }


                all_probs(r,s) = sum_prob;

            }
            if(all_em_stats[r].stat != all_stats_same[r] || all_em_stats[r].stat_diff!=all_stats_diff[r]){
                cout << "DIFF!!! outer " << endl;
            }
        }

        MaximizationStep();


        for (size_t r = 0; r < rate_count; ++r) {
            auto em_stat = all_em_stats[r];
            if(em_stat.stat != all_stats_same[r] || em_stat.stat_diff!=all_stats_diff[r]){
                cout << "DIFF!!! max" << em_stat.stat << "\t" <<all_stats_same[r] << "\t" << em_stat.stat_diff << "\t" << all_stats_diff[r] <<
                        endl;

            }
            double sum_stat = all_stats_diff[r] + all_stats_same[r];
            double new_exp_beta = all_stats_diff[r] / sum_stat;
            double new_mu = mutation_prob.ConvertExpBetaToMu(new_exp_beta);

            double new_one_minus_exp_beta = all_stats_same[r] / sum_stat;
            double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta-1));
//            cout.setprecision(10);

//            proportion[r] = all_probs[r]/(all_probs[0]+all_probs[1]);
            printf("======= NEM_MU_r: %zu \tMu: %.5e %.5f =expBeta=  %.5f %.5f \t =prop= %.5f %.5f "
                            "\t =stat= %.5f %.5f \n" ,r,
                    new_mu ,new_one_minus_mu,
                    new_exp_beta, new_one_minus_exp_beta,
                    proportion[0], proportion[1],
                    all_stats_same[0], all_stats_same[1]);

            parameters[r] = new_mu;

        }


    }

}


void EM::InitialiseParameters() {
    double lower_bound = 1e-100;
    double upper_bound = 1;
    parameters = vector<double>(num_category);
    if (num_category == 2){
        parameters = {upper_bound, lower_bound};
    }
    else{
        cout << "Not yet implemented for more than 2 categories" << endl;
        exit(222);
        //TODO: Should throw exception instead of exit, this will do for now
    }

}

void EM::InitialiseProportion() {
    double default_proportion = 1.0/num_category;
    proportion = vector<double> (num_category, default_proportion);
}

void EM::CalculateProportion() {

    auto sum_col = all_probs.colwise().sum();
    auto prop_col = all_probs.row(0) / sum_col;
    double proportion_sum = prop_col.sum();
    proportion[0] = proportion_sum / site_count;
    proportion[1] = 1 - proportion[0];

    if (num_category != 2) {
        cout << "Not yet implemented for more than 2 categories" << endl;
        exit(222);
    }

}

void EM::MaximizationStep() {
    CalculateProportion();


}
