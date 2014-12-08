/*
 * em_algorithm_mutation.h
 *
 *  Created on: 12/3/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_ALGORITHM_MUTATION_H_
#define EM_ALGORITHM_MUTATION_H_


#include <bits/unique_ptr.h>
#include <memory>
//#include <stddef.h>
//#include "site_prob.h"

#include "em_algorithm.h"

#include "em_model_mutation.h"
#include "em_data_mutation.h"
#include "em_summary_stat_mutation.h"


class EmAlgorithmMutation : public EmAlgorithm{

public:

//    EM(vector<EmData> em_data, EmModel &em_model);

//    EM(int num_category, vector<SiteProb> em_data0, EvolutionModel &em_model0, EmModel &m);


//    EM(vector<EmDataMutation> site_data, EmModelMutation &evo_model);

//    EM(int num_category0, vector<EmData> em_data0, EmModel &em_model0);

//    EM(int num_category0, vector<SiteProb> em_data0, EmModel &em_model0);

//    EM(int num_category0, vector<SiteProb> em_data0);

//    EM(int num_category0, vector<SiteProb> em_data0, EvolutionModel &em_model0, EmModelMutation &m);

//    EM(int num_category0, vector<SiteProb> &em_data0, EvolutionModel &em_model0, vector<EmData*> &d, EmModel &m);

    EmAlgorithmMutation(int num_category0, vector<SiteProb> &em_data0, EvolutionModel &em_model0, vector<unique_ptr<EmData>> &d_ptr, EmModel &m);

    EmAlgorithmMutation(vector<unique_ptr<EmData>> &d_ptr, vector<unique_ptr<EmModel>> &m);

    virtual ~EmAlgorithmMutation();

    void Run();






private:


protected:

    virtual void InitialiseParameters();
    virtual void InitialiseSummaryStat();





private:
    //TODO These will be gone soon
    std::vector<double> all_stats_same;
    std::vector<double> all_stats_diff;

    EvolutionModel *em_model_old;

    vector<double> parameters_old;

    Eigen::ArrayXXd all_probs_old;

    vector<SiteProb> em_data_old;
    void oldEStep();

    void oldMStep();

};


#endif //EM_ALGORITHM_MUTATION_H_

