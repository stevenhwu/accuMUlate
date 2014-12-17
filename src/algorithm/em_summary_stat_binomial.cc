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


//    printf("===================================== NEM_Theta: %.5f %.5f %.5f\n",
//            stat[0], new_theta, sum_stat);


    return new_theta;
}
/*

void EmSummaryStatBinomial::biTest(double a){

}

*/