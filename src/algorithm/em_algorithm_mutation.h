/*
 * em_algorithm_mutation.h
 *
 *  Created on: 12/3/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_ALGORITHM_MUTATION_H_
#define EM_ALGORITHM_MUTATION_H_


//#include <bits/unique_ptr.h>
#include <memory>
//#include <stddef.h>
//#include "site_prob.h"

#include "em_algorithm.h"

#include "em_model_mutation.h"



class EmAlgorithmMutation : public EmAlgorithm{

public:

    EmAlgorithmMutation(int num_category0,  EmModelMutation &em_model);

    EmAlgorithmMutation(std::vector<std::unique_ptr<EmModel>> &model_ptr);
    virtual ~EmAlgorithmMutation();

    void Run();

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

