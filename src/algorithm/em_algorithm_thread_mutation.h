/*
 * em_algorithm_mutation.h
 *
 *  Created on: 12/3/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_ALGORITHM_THREAD_MUTATION_H_
#define EM_ALGORITHM_THREAD_MUTATION_H_

#include "em_algorithm_thread.h"


class EmAlgorithmThreadMutation : public EmAlgorithmMultiThreading{

public:


    EmAlgorithmThreadMutation(MutationModelMultiCategories &model_multi, uint32_t thread_count);

    virtual ~EmAlgorithmThreadMutation();

    void RunEM();



protected:


//    virtual void InitialiseParameters();
    virtual void InitialiseSummaryStat();

    virtual void ExpectationStepCustom(size_t data_index, size_t category_index, double &sum_prob, std::vector<double> &temp_stat);

};


#endif //EM_ALGORITHM_THREAD_MUTATION_H_

