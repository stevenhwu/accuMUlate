/*
 * em_algorithm_mutation.h
 *
 *  Created on: 12/3/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_ALGORITHM_MUTATION_H_
#define EM_ALGORITHM_MUTATION_H_


#include <memory>
//#include <stddef.h>
//#include "site_prob.h"

#include "em_algorithm.h"


class EmAlgorithmMutation : public EmAlgorithm{

public:

    EmAlgorithmMutation(int num_category0,  EmModelMutation &em_model);

    EmAlgorithmMutation(std::vector<std::unique_ptr<EmModel>> &model_ptr);

    virtual ~EmAlgorithmMutation();

    void RunEM();


private:


protected:

//    virtual void InitialiseParameters();
    virtual void InitialiseSummaryStat();

    virtual void ExpectationStepCustom(size_t data_index, size_t category_index, double &sum_prob, std::vector<double> &temp_stat);

};


#endif //EM_ALGORITHM_MUTATION_H_

