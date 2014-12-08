/*
 * em_summary_stat_mutation.cc
 *
 *  Created on: 12/8/14
 *      Author: Steven Wu
 */


#include "em_summary_stat_mutation.h"


double EmSummaryStatMutation::MaximiseStats() {

    double sum_stat = stat[0] + stat[1];
    double new_exp_beta = stat[1] / sum_stat;
//    double new_mu = mutation_prob.ConvertExpBetaToMu(new_exp_beta);

    double new_one_minus_exp_beta = stat[0] / sum_stat;
//    double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta - 1));


    printf("======= NEM_MU_r:\t\t\t\t\t\t \t\t\t=expBeta=  %.5f %.5f \n",
            new_exp_beta, new_one_minus_exp_beta
    );
    return new_exp_beta;
    /*
    double sum_stat = all_stats_diff[r] + all_stats_same[r];
    double new_exp_beta = all_stats_diff[r] / sum_stat;
    double new_mu = mutation_prob.ConvertExpBetaToMu(new_exp_beta);

    double new_one_minus_exp_beta = all_stats_same[r] / sum_stat;
    double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta - 1));
//            cout.setprecision(10);

//            proportion[r] = all_probs[r]/(all_probs[0]+all_probs[1]);


    all_em_stats[r]->MaximiseStats();



    printf("======= NEM_MU_r: %zu \tMu: %.5e %.5f =expBeta=  %.5f %.5f \t =prop= %.5f %.5f "
                    "\t =stat= %.5f %.5f \n", r,
            new_mu, new_one_minus_mu,
            new_exp_beta, new_one_minus_exp_beta,
            proportion[0], proportion[1],
            all_stats_same[0], all_stats_same[1]);
    */
}
