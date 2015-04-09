/*
 * em_algorithm_mutation.h
 *
 *  Created on: 12/3/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_ALGORITHM_MUTATION_V1_H_
#define EM_ALGORITHM_MUTATION_V1_H_



#include <memory>
//#include <stddef.h>
//#include "site_prob.h"

#include "em_algorithm.h"

#include "em_model_mutation_v1.h"
#include "em_data_mutation_v1.h"
#include "em_summary_stat_mutation.h"


class EmAlgorithmMutationV1 : public EmAlgorithm{

public:

    EmAlgorithmMutationV1(std::vector<std::unique_ptr<EmData>> &data_ptr,
            std::vector<std::unique_ptr<EmModel>> &em_model_ptr);

    EmAlgorithmMutationV1(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr, EmModelMutationV1 &em_model);

    virtual ~EmAlgorithmMutationV1();

    void RunEM();

    void Run2();

private:


protected:

    virtual void InitialiseParameters();
    virtual void InitialiseSummaryStat();


    virtual void ExpectationStepCustom(size_t data_index, size_t category_index, double &sum_prob, std::vector<double> &temp_stat);



};


#endif //EM_ALGORITHM_MUTATION_H_

