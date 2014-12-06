/*
 * em_algorithm.h
 *
 *  Created on: 12/3/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_ALGORITHM_H_
#define EM_ALGORITHM_H_


#include <bits/unique_ptr.h>
#include <stddef.h>
#include "site_prob.h"

#include "em_model_mutation.h"
#include "em_data_mutation.h"
#include "em_summary_stat_mutation.h"


class EM {

public:

//    EM(vector<EmData> em_data, EmModel &em_model);

//    EM(int num_category, vector<SiteProb> em_data0, EvolutionModel &em_model0, EmModel &m);


//    EM(vector<EmDataMutation> site_data, EmModelMutation &evo_model);

//    EM(int num_category0, vector<EmData> em_data0, EmModel &em_model0);

//    EM(int num_category0, vector<SiteProb> em_data0, EmModel &em_model0);

//    EM(int num_category0, vector<SiteProb> em_data0);

//    EM(int num_category0, vector<SiteProb> em_data0, EvolutionModel &em_model0, EmModelMutation &m);

//    EM(int num_category0, vector<SiteProb> &em_data0, EvolutionModel &em_model0, vector<EmData*> &d, EmModel &m);

    EM(int num_category0, vector<SiteProb> &em_data0, EvolutionModel &em_model0, vector<unique_ptr<EmData>> &d_ptr, EmModel &m);

    virtual ~EM();

    void Run();
    void RunOld();
protected:
    void ExpectationStep ();

    void MaximizationStep();


private:
    //TODO These will be gone soon
    vector<double> all_stats_same;
    vector<double> all_stats_diff;


private:

    int num_category;

    vector<SiteProb> em_data_old;
    EvolutionModel *em_model_old;

    vector<double> parameters_old;
//    vector<EmData*> em_data;
    vector<unique_ptr<EmData>> *em_data_ptr;
    EmModel *em_model;

    vector<double> parameters;
    vector<double> proportion;


    size_t rate_count;
    size_t site_count;
    size_t em_count;

    void InitialiseParameters();

    void Init();

    void InitialiseProportion();

    Eigen::ArrayXXd all_probs;
    Eigen::ArrayXXd all_probs_test;


    void CalculateProportion();

    vector<unique_ptr<EmSummaryStat>> all_em_stats;


    void oldEStep();

//    int parameters_old;
    void oldMStep(MutationProb &mutation_prob);
};


#endif //EM_ALGORITHM_H_

