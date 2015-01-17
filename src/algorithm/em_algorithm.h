/*
 * em_algorithm.h
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */


#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <stddef.h>
#include <Eigen/Dense>
#include "em_model.h"
#include "em_data.h"

#ifndef EM_ALGORITHM_H_
#define EM_ALGORITHM_H_

extern const double EM_CONVERGE_THRESHOLD;

class EmAlgorithm {


public:

    EmAlgorithm(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr, EmModel &em_model0);

    EmAlgorithm(std::vector<std::unique_ptr<EmData>> &data_ptr, std::vector<std::unique_ptr<EmModel>> &model_ptr);

//    EmAlgorithm() ;

//    EmAlgorithm(int category_count);

    EmAlgorithm(std::vector<std::unique_ptr<EmModel>> &model_ptr);

    virtual ~EmAlgorithm() {
    }


    virtual void Run() = 0;

    virtual std::vector<double> GetProportion();

    virtual std::vector<double> GetParameters();

    void PrintSummary();

protected:


    size_t num_category;
    size_t site_count;
    size_t max_ite_count;

    std::vector<std::unique_ptr<EmData>> *em_data_ptr;
    std::vector<std::unique_ptr<EmModel>> *em_model_ptr;

    std::vector<std::unique_ptr<EmModel>> em_model;
    EmModel *em_model0;

    std::vector<std::unique_ptr<EmSummaryStat>> all_em_stats;

    std::vector<std::unique_ptr<EmSummaryStat>> em_stat_local_ptr;
    std::unique_ptr<EmSummaryStat> em_stat_local_single;


    std::vector<double> parameters;
    std::vector<double> proportion;
    Eigen::ArrayXXd all_probs;

    std::vector<std::vector<double>> temp_stats;
    std::vector<double> cache_parameters;


    void InitWithData();

    void InitWithModel();

    void ExpectationStepModel();

    void ExpectationStepModelPtr();

    void MaximizationStep();

    void CalculateProportion();

    virtual void InitialiseProportion();

    virtual void InitialiseParameters() = 0;

    virtual void InitialiseSummaryStat() = 0;

    virtual void ExpectationStepCustom(size_t data_index, size_t category_index,
            double &sum_prob, std::vector<double> &temp_stat) = 0;

    bool EmStoppingCriteria(int ite);


};


#endif //EM_ALGORITHM_H_
