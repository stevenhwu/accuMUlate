/*
 * em_algorithm.h
 *
 *  Created on: 12/3/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_ALGORITHM_H_
#define EM_ALGORITHM_H_


#include "site_prob.h"

#include "em_model_mutation.h"
#include "em_data_mutation.h"
#include "em_summary_stat_mutation.h"


class EM {

public:

    EM(vector<EmData> em_data, EmModel &em_model);

    EM(int num_category, vector<SiteProb> em_data0, EvolutionModel &em_model0);


    EM(vector<EmDataMutation> site_data, EmModelMutation &evo_model);

    virtual ~EM();

    void Run();

protected:
    void ExpectationStep ();

    void MaximizationStep();


private:
    //TODO These will be gone soon
    vector<double> all_stats_same;
    vector<double> all_stats_diff;


private:

    int num_category;

    vector<SiteProb> em_data;
    EvolutionModel *em_model;

//    vector<EmData> em_data;
//    EmModel *em_model;

    vector<double> parameters;
    vector<double> proportion;


    size_t rate_count;
    size_t site_count;
    size_t em_count;

    void InitialiseParameters();

    void Init();

    void InitialiseProportion();

    Eigen::ArrayXXd all_probs;


    void CalculateProportion();

    vector<EmSummaryStatMutation> all_em_stats;
};


#endif //EM_ALGORITHM_H_

