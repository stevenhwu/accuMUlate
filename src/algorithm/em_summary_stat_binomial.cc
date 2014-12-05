#include <ostream>
#include <iostream>
#include "em_summary_stat_binomial.h"



void EmSummaryStatBinomial::print() {
//    std::cout << "EMSumStat2: " << stat << "\t" << stat_diff << std::endl;
//    EmSummaryStat::stat;
}

EmSummaryStatBinomial::EmSummaryStatBinomial() {

}

void EmSummaryStatBinomial::SetStats(double stat_same0, double stat_diff0) {
    stat = stat_same0;
    stat_diff = stat_diff0;
}

void EmSummaryStatBinomial::Reset() {
//    stat = 0;
    EmSummaryStat::Reset();
    stat_diff = 0;


}


void EmSummaryStatBinomial::UpdateSumWithProportion(double &d, EmSummaryStatBinomial mutation) {
    stat += d* mutation.stat;
    stat_diff += d* mutation.stat_diff;


}

void EmSummaryStatBinomial::SetStats(std::vector<double> stats) {

    stat = stats[0];
    stat_diff = stats[1];
    std::cout << "SetStat in EMSumStat_MUTATINO!!: " << stat  << "\t" << stat_diff << std::endl;
}


void EmSummaryStatBinomial::biTest(double a){

}