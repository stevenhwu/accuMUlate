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

//    double new_one_minus_exp_beta = stat[0] / sum_stat;
//    double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta - 1));


//    printf("======= NEM_MU_r:\t\t\t\t\t\t \t\t\t=expBeta=  %.5f %.5f \n",
//            new_exp_beta, new_one_minus_exp_beta
//    );
    return new_exp_beta;

}
