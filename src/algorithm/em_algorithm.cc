#include "em_algorithm.h"
#include "em_data_mutation.h"
#include "em_model_mutation.h"

#include <memory>
#include <stddef.h>
//EM::EM(int num_category0, vector<SiteProb> em_data0, EvolutionModel &em_model0) : num_category(num_category0), em_data_old(em_data0), em_model_old(&em_model0){
//
//    cout << "EM Constructor\n";
//
//    if (num_category != 2) {
//        cout << "Not yet implemented for more than 2 categories" << endl;
//        exit(222);
//    }
//
//    site_count = em_data_old.size();
//    rate_count = num_category;
//    em_count = 10;
//
//    Init();
//
//
//
//}

//EM::EM(int num_category0, vector<SiteProb> em_data0, EmModel &em_model0) : num_category(num_category0), em_data_old(em_data0) {//}, em_model(&em_model0){

//std::vector<std::unique_ptr<EmData>>
EM::EM(int num_category0, vector<SiteProb> &em_data0, EvolutionModel &em_model0, vector<std::unique_ptr<EmData>> &d_ptr, EmModel &m)
        : num_category(num_category0), em_data_old(em_data0), em_model_old(&em_model0), em_data_ptr(&d_ptr), em_model(&m) {

//    cout << "\n============\nEM  half way Constructor\n";

//    em_model->UpdateParameter(2);
    if (num_category != 2) {
        cout << "Not yet implemented for more than 2 categories" << endl;
        exit(222);
    }

    site_count = em_data_ptr->size();
    rate_count = num_category;
    em_count = 100;

//    auto site2 = em_data[0];
//    EmData *site3 = em_data[0];
//
//    site2->Test(4);
//    site3->Test(5);
    double lower_bound = 1e-12;
    double upper_bound = 1e-1;
    all_probs_test = Eigen::ArrayXXd::Zero(num_category, site_count);
    parameters_old = vector<double>(num_category);
    if (num_category == 2) {
        parameters_old = {upper_bound, lower_bound};
    }
    Init();
    cout << "===========\nDone Constructor smart pointer\n";

}


//
//EM::EM(int num_category0, vector<SiteProb> &em_data0, EvolutionModel &em_model0, vector<EmData *> &d, EmModel &m)
//        : num_category(num_category0), em_data_old(em_data0), em_model_old(&em_model0), em_data(d), em_model(&m) {
//
////    cout << "\n============\nEM  half way Constructor\n";
//
////    em_model->UpdateParameter(2);
//    if (num_category != 2) {
//        cout << "Not yet implemented for more than 2 categories" << endl;
//        exit(222);
//    }
//
//    site_count = em_data.size();
//    rate_count = num_category;
//    em_count = 3;
//
////    auto site2 = em_data[0];
////    EmData *site3 = em_data[0];
////
////    site2->Test(4);
////    site3->Test(5);
//
//
//    Init();
//    cout << "===========\nDone Constructor raw pointer\n";
//}


//EM::EM(int num_category0, vector<EmData> em_data0, EmModel &em_model0) : num_category(num_category0), em_data(em_data0), em_model(&em_model0){
//
//    cout << "EM Constructor\n";
//    if (num_category != 2) {
//        cout << "Not yet implemented for more than 2 categories" << endl;
//        exit(222);
//    }
//
//    site_count = em_data_old.size();
//    rate_count = num_category;
//    em_count = 10;
//
//    Init();
//}


EM::~EM() {

}

void EM::Init() {

    InitialiseParameters();
    InitialiseProportion();
    all_probs = Eigen::ArrayXXd::Zero(num_category, site_count);
    all_stats_same = vector<double>(num_category, 0);
    all_stats_diff = vector<double>(num_category, 0);

//    all_em_stats = vector<EmSummaryStatMutation>(rate_count);
    for (int i = 0; i < rate_count; ++i) {
        all_em_stats.emplace_back(new EmSummaryStatMutation());
    }

}


void EM::ExpectationStep() {

    EmSummaryStat em_stat_local = EmSummaryStatMutation();

    for (size_t r = 0; r < rate_count; ++r) {
        em_model->UpdateParameter(parameters[r]);

        all_em_stats[r]->Reset();

        for (size_t s = 0; s < site_count; ++s) {
            cout << "em_count: "  << " rate_count: " << r << " site: " << s << endl;
            double sum_prob = 0;

            (*em_data_ptr)[s]->UpdateEmModel(em_model);

            (*em_data_ptr)[s]->UpdateSummaryStat(sum_prob, em_stat_local);

            all_em_stats[r]->UpdateSumWithProportion(proportion[r], em_stat_local);

            all_probs(r, s) = sum_prob;

        }
    }
}

void EM::Run() {

    MutationProb mutation_prob = em_model_old->GetMutationProb();//FIXME: What to do here?? maybe in EmSummayrStat??

    EmSummaryStatMutation em_stat_local;

    for (size_t i = 0; i < em_count; ++i) {
        cout << "EM ite: " << i << endl;
        ExpectationStep();

        oldEStep();
//        for (int r = 0; r < rate_count; ++r) {
//            if (all_em_stats[r]->GetStat(0) != all_stats_same[r] || all_em_stats[r]->GetStat(1) != all_stats_diff[r]) {
//                cout << "DIFF!!! inner " << all_em_stats[r]->GetStat(0) << "\t" << all_stats_same[r] << "\t" << all_em_stats[r]->GetStat(1) << "\t" << all_stats_diff[r] << endl;
//                exit(99);
//            }
//        }
//        Eigen::ArrayXXd complex1 = (all_probs - all_probs_test);
//        double sumaa = complex1.sum();
//        if (sumaa != 0) {
//            cout << all_probs;
//            cout << "=============================\n";
//            cout << all_probs_test;
//            cout << sumaa << endl;
//            exit(199);
//        }

        MaximizationStep();


        for (size_t r = 0; r < rate_count; ++r) {
//            parameters[r] = new_mu;
            parameters[r] = mutation_prob.ConvertExpBetaToMu(parameters[r]);

        }
        oldMStep(mutation_prob);
        if(parameters[0] != parameters_old[0] || parameters[1] != parameters_old[1]){
            cout << "DIFF param: " << parameters[0] << "\t" <<  parameters_old[0]<< "\t" <<  parameters[1] << "\t" <<  parameters_old[1] << endl;
            exit(144);
        }

    }

}

void EM::oldMStep(MutationProb &mutation_prob) {
    for (size_t r = 0; r < rate_count; ++r) {

//            all_em_stats[r]->MaximiseStats();
            double sum_stat = all_stats_diff[r] + all_stats_same[r];
            double new_exp_beta = all_stats_diff[r] / sum_stat;
            double new_mu = mutation_prob.ConvertExpBetaToMu(new_exp_beta);

            double new_one_minus_exp_beta = all_stats_same[r] / sum_stat;
            double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta - 1));
//            cout.setprecision(10);

//            proportion[r] = all_probs[r]/(all_probs[0]+all_probs[1]);


            printf("======= NEM_MU_r: %zu \tMu: %.5e %.5f =expBeta=  %.5f %.5f \t =prop= %.5f %.5f "
                            "\t =stat= %.5f %.5f \n", r,
                    new_mu, new_one_minus_mu,
                    new_exp_beta, new_one_minus_exp_beta,
                    proportion[0], proportion[1],
                    all_stats_same[0], all_stats_same[1]);

            parameters_old[r] = new_mu;

        }
}

void EM::oldEStep() {
    for (size_t r = 0; r < rate_count; ++r) {

            em_model_old->UpdateMu(parameters[r]);

            all_stats_same[r] = 0;
            all_stats_diff[r] = 0;

//            all_em_stats[r]->Reset();
//            em_model->UpdateParameter(parameters[r]+0.05);

            for (size_t s = 0; s < site_count; ++s) {
                cout << "em_E: "  << " rate_count: " << r << " site: " << s << endl;
                double sum_prob = 0;
                auto site_old = em_data_old[s];
                site_old.UpdateModel(*em_model_old);
                double stat_same = 0;
                double stat_diff = 0;
                site_old.CalculateAncestorToDescendant(sum_prob, stat_same, stat_diff);
                all_stats_same[r] += proportion[r] * stat_same;
                all_stats_diff[r] += proportion[r] * stat_diff;


                all_probs_test(r, s) = sum_prob;

            }
//            if(all_em_stats[r].stat != all_stats_same[r] || all_em_stats[r].stat_diff!=all_stats_diff[r]){
//                cout << "DIFF!!! outer " << endl;
//            }
        }
}


void EM::RunOld() {

    MutationProb mutation_prob = em_model_old->GetMutationProb();


//    for (size_t s = 0; s < site_count; ++s) {
//        EmSummaryStatMutation e;
////        e.stat = 1;
////        e.stat_diff = 2;
//        all_em_stats.push_back(e);
//    }



    for (size_t i = 0; i < em_count; ++i) {
        cout << "EM ite: " << i << endl;

        for (size_t r = 0; r < rate_count; ++r) {
            em_model_old->UpdateMu(parameters[r]);
//            mutation_prob.UpdateMu(parameters[r]);

            all_stats_same[r] = 0;
            all_stats_diff[r] = 0;
            //replace with
//            auto em_stat = all_em_stats[r];
            all_em_stats[r]->Reset();

//            em_stat.print();
//            all_probs[r][s] = 0;
//          Array10D site_stat[site_count];
            for (size_t s = 0; s < site_count; ++s) {
                double sum_prob = 0;

                auto site = em_data_old[s];

                site.UpdateModel(*em_model_old);

                double stat_same = 0;
                double stat_diff = 0;
                //relpace with
                EmSummaryStatMutation local;


                site.CalculateAncestorToDescendant(sum_prob, stat_same, stat_diff);


                all_stats_same[r] += proportion[r] * stat_same;
                all_stats_diff[r] += proportion[r] * stat_diff;
                //replace with
                local.SetStats(vector<double>{ stat_same, stat_diff});
                all_em_stats[r]->UpdateSumWithProportion(proportion[r], local);


//                if(all_em_stats[r].stat != all_stats_same[r] || all_em_stats[r].stat_diff!=all_stats_diff[r]){
//                    cout << "DIFF!!! inner " << endl;
//                }


                all_probs(r, s) = sum_prob;

            }
//            if(all_em_stats[r].stat != all_stats_same[r] || all_em_stats[r].stat_diff!=all_stats_diff[r]){
//                cout << "DIFF!!! outer " << endl;
//            }
        }

        MaximizationStep();

        for (size_t r = 0; r < rate_count; ++r) {
//            parameters[r] = new_mu;
            parameters[r] = mutation_prob.ConvertExpBetaToMu(parameters[r]);

        }
        for (size_t r = 0; r < rate_count; ++r) {
//            auto em_stat = all_em_stats[r];
//            if(em_stat.stat != all_stats_same[r] || em_stat.stat_diff!=all_stats_diff[r]){
//                cout << "DIFF!!! max" << em_stat.stat << "\t" <<all_stats_same[r] << "\t" << em_stat.stat_diff << "\t" << all_stats_diff[r] <<
//                        endl;
//
//            }
            double sum_stat = all_stats_diff[r] + all_stats_same[r];
            double new_exp_beta = all_stats_diff[r] / sum_stat;
            double new_mu = mutation_prob.ConvertExpBetaToMu(new_exp_beta);

            double new_one_minus_exp_beta = all_stats_same[r] / sum_stat;
            double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta - 1));
//            cout.setprecision(10);

//            proportion[r] = all_probs[r]/(all_probs[0]+all_probs[1]);
            printf("======= NEM_MU_r: %zu \tMu: %.5e %.5f =expBeta=  %.5f %.5f \t =prop= %.5f %.5f "
                            "\t =stat= %.5f %.5f \n", r,
                    new_mu, new_one_minus_mu,
                    new_exp_beta, new_one_minus_exp_beta,
                    proportion[0], proportion[1],
                    all_stats_same[0], all_stats_same[1]);

            parameters_old[r] = new_mu;

        }


    }

}


void EM::InitialiseParameters() {
    double lower_bound = 1e-12;
    double upper_bound = 1e-1;
    parameters = vector<double>(num_category);
    if (num_category == 2) {
        parameters = {upper_bound, lower_bound};
    }
    else {
        cout << "Not yet implemented for more than 2 categories" << endl;
        exit(222);
        //TODO: Should throw exception instead of exit, this will do for now
    }

}

void EM::InitialiseProportion() {
    double default_proportion = 1.0 / num_category;
    proportion = vector<double>(num_category, default_proportion);
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


    for (size_t r = 0; r < rate_count; ++r) {

        double new_parameter = all_em_stats[r]->MaximiseStats();
//        double new_mu = mutation_prob.ConvertExpBetaToMu(new_exp_beta);

//        double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta - 1));


        parameters[r] = new_parameter;

    }


//    printf("======= NEM_MU_r: %zu \tMu: %.5e %.5f =expBeta=  %.5f %.5f \t =prop= %.5f %.5f "
//                    "\t =stat= %.5f %.5f \n", r,
//            new_mu, new_one_minus_mu,
//            new_exp_beta, new_one_minus_exp_beta,
//            proportion[0], proportion[1],
//            all_stats_same[0], all_stats_same[1]

}

