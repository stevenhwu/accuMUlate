
#include "em_summary_stat_binomial.h"



EmSummaryStatBinomial::EmSummaryStatBinomial() : EmSummaryStat(EM_SUMMARY_STAT_BINOMIAL_STATS_COUNT) {

}


double EmSummaryStatBinomial::MaximiseStats() {

    double sum_stat = stat[0] + stat[1];
    double new_theta = stat[1] / sum_stat;
    return new_theta;
}
