
#include <iostream>
#include "em_summary_stat_binomial.h"


//
//void EmSummaryStatBinomial::print() {
////    std::cout << "EMSumStat2: " << stat << "\t" << stat_diff << std::endl;
////    EmSummaryStat::stat;
//}

EmSummaryStatBinomial::EmSummaryStatBinomial() {
    stat_count = 2;
    stat = std::vector<double> (stat_count);

}


void EmSummaryStatBinomial::biTest(double a){

}