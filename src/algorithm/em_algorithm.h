/*
 * em_algorithm.h
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */


#pragma once

#include <memory>
#include <vector>
#include <stddef.h>
#include <Eigen/Dense>
#include "em_model.h"
#include "em_data.h"

#ifndef EM_ALGORITHM_H_
#define EM_ALGORITHM_H_


class EmAlgorithm {


    void CalculateProportion();

public:

    EmAlgorithm(int num_category0, std::vector <std::unique_ptr<EmData>> &data_ptr, EmModel &em_model0);

//    EmAlgorithm() ;

//    EmAlgorithm(int category_count);

    virtual ~EmAlgorithm() {
    }

    virtual std::vector<double> GetParameters();


    virtual void Run() = 0;

protected:


    size_t num_category;
    size_t site_count;
    size_t em_count;

    std::vector<std::unique_ptr<EmData>> *em_data_ptr;
    std::vector<std::unique_ptr<EmModel>> *em_model_ptr;
    EmModel *em_model;

    std::vector<std::unique_ptr<EmSummaryStat>> all_em_stats;
    std::unique_ptr<EmSummaryStat> em_stat_local;

    std::vector<double> parameters;
    std::vector<double> proportion;
    Eigen::ArrayXXd all_probs;




    void Init();


    void ExpectationStep ();

    void MaximizationStep();

    virtual void InitialiseProportion();



    virtual void InitialiseParameters() = 0;

    virtual void InitialiseSummaryStat() = 0;

    void ExpectationStep2();

    EmAlgorithm(std::vector<std::unique_ptr<EmData>> &data_ptr, std::vector<std::unique_ptr<EmModel>> &model_ptr);
};


#endif //EM_ALGORITHM_H_
