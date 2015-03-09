/*
 * em_algorithm_mutation.h
 *
 *  Created on: 12/3/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_ALGORITHM_MUTATION_V1_H_
#define EM_ALGORITHM_MUTATION_V1_H_


//#include <bits/unique_ptr.h>
#include <memory>
//#include <stddef.h>
//#include "site_prob.h"

#include "em_algorithm.h"

#include "em_model_mutation_v1.h"
#include "em_data_mutation_v1.h"
#include "em_summary_stat_mutation_v1.h"


class EmAlgorithmMutationV1 : public EmAlgorithm{

public:

//    EM(vector<EmData> em_data, EmModel &em_model);

//    EM(int num_category, vector<SiteProb> em_data0, EvolutionModel &em_model0, EmModel &m);

//    EM(vector<EmDataMutationV1> site_data, EmModelMutationV1 &evo_model);

//    EM(int num_category0, vector<EmData> em_data0, EmModel &em_model0);

//    EM(int num_category0, vector<SiteProb> em_data0, EmModel &em_model0);

//    EM(int num_category0, vector<SiteProb> em_data0);

//    EM(int num_category0, vector<SiteProb> em_data0, EvolutionModel &em_model0, EmModelMutationV1 &m);

//    EM(int num_category0, vector<SiteProb> &em_data0, EvolutionModel &em_model0, vector<EmData*> &d, EmModel &m);

    EmAlgorithmMutationV1(int num_category0, std::vector<SiteProb> &em_data0, EvolutionModel &em_model0, std::vector<std::unique_ptr<EmData>> &data_ptr, EmModel &em_model);

    EmAlgorithmMutationV1(std::vector<std::unique_ptr<EmData>> &data_ptr, std::vector<std::unique_ptr<EmModel>> &em_model_ptr);

    EmAlgorithmMutationV1(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr, EmModelMutationV1 &em_model);

    virtual ~EmAlgorithmMutationV1();

    void Run();


    void Run2();

private:


protected:

    virtual void InitialiseParameters();
    virtual void InitialiseSummaryStat();


    virtual void ExpectationStepCustom(size_t data_index, size_t category_index, double &sum_prob, std::vector<double> &temp_stat);


private:
    //TODO These will be gone soon
    std::vector<double> all_stats_same;
    std::vector<double> all_stats_diff;

    std::vector<SiteProb> em_data_old;
    EvolutionModel *em_model_old;

    std::vector<double> parameters_old;
    Eigen::ArrayXXd all_probs_old;



};


#endif //EM_ALGORITHM_MUTATION_H_

