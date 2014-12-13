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

    EmAlgorithmMutation(int num_category0, std::vector<SiteProb> &em_data0, EvolutionModel &em_model0, std::vector<std::unique_ptr<EmData>> &data_ptr, EmModel &em_model);

    EmAlgorithmMutation(std::vector<std::unique_ptr<EmData>> &data_ptr, std::vector<std::unique_ptr<EmModel>> &em_model_ptr);

    EmAlgorithmMutation(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr, EmModelMutation &em_model);

    virtual ~EmAlgorithmMutation();

    void Run();


    void Run2();

private:


protected:

    virtual void InitialiseParameters();
    virtual void InitialiseSummaryStat();





private:
    //TODO These will be gone soon
    std::vector<double> all_stats_same;
    std::vector<double> all_stats_diff;

    std::vector<SiteProb> em_data_old;
    EvolutionModel *em_model_old;
    std::vector<double> parameters_old;

    Eigen::ArrayXXd all_probs_old;


    void oldEStep();

    void oldMStep();


};


#endif //EM_ALGORITHM_MUTATION_H_

