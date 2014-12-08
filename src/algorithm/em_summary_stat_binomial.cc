#include "em_summary_stat_binomial.h"


//
//void EmSummaryStatBinomial::print() {
////    std::cout << "EMSumStat2: " << stat << "\t" << stat_diff << std::endl;
////    EmSummaryStat::stat;
//}

EmSummaryStatBinomial::EmSummaryStatBinomial() : EmSummaryStat(EM_SUMMARY_STAT_BINOMIAL_STATS_COUNT) {

}


double EmSummaryStatBinomial::MaximiseStats() {

    double sum_stat = stat[0] + stat[1];
    double new_theta = stat[0] / sum_stat;
//    double new_mu = mutation_prob.ConvertExpBetaToMu(new_exp_beta);

    double new_one_minus_theta = stat[1] / sum_stat;
//    double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta - 1));


    printf("===================================== NEM_Theta: %.5f %.5f \n",
            new_theta, new_one_minus_theta
    );
    return new_theta;
}
/*

void EmSummaryStatBinomial::biTest(double a){

}

*/