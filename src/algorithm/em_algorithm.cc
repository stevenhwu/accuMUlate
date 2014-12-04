#include "em_algorithm.h"

EM::EM(vector<SiteProb> site_data, EvolutionModel &evo_model): site_data(site_data), evo_model(&evo_model){

    site_count = site_data.size();
}



EM::~EM() {

}

void EM::Single() {

//    int all_prob = 0;

    const size_t cat = 2;

    double muArray[cat];
    muArray[0] = 1e-10;
    muArray[1] = 1;


    vector<double> proportion;// {0.5,0.5};
    vector<double> all_stats_same;//[cat];
    vector<double> all_stats_diff;//[cat];

    MutationProb mutation_prob = evo_model-> GetMutationProb();
//    size_t em_count = 50;
//    size_t rate_count = 2; //2;
    double all_prob[cat][site_count];



    for (size_t i = 0; i < em_count; ++i) {



        for (size_t r = 0; r < rate_count; ++r) {
            evo_model->UpdateMu(muArray[r]);
//            mutation_prob.UpdateMu(muArray[r]);

            all_stats_same[r] = 0;
            all_stats_diff[r] = 0;
//            all_prob[r][s] = 0;
//          Array10D site_stat[site_count];
            for (size_t s = 0; s < site_count; ++s) {

                auto site = site_data[s];
//                t.UpdateTransitionMatrix(model);
//                t.UpdateMuProb(mutation_prob);
                site.UpdateModel(*evo_model);
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
        for (size_t s = 0; s < site_count; ++s) {
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

void EMData2::UpdateLikelihood(double prob, EMSummaryStat &summaryStat) {
    cout << "call EMDATA2" << endl;
    summaryStat.print();
}

void EMSummaryStat2::print() {
    cout << "EMSumStat2: " << stat1 << "\t" << stat2 <<endl;
}
