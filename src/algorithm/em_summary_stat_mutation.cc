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
//    double new_exp_beta = (0.0+stat[1]) / (0.0+sum_stat);
    return new_exp_beta;

}
