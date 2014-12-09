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



public:

    EmAlgorithm(int num_category0, std::vector <std::unique_ptr<EmData>> &data_ptr, EmModel &em_model0);

    EmAlgorithm(std::vector<std::unique_ptr<EmData>> &data_ptr, std::vector<std::unique_ptr<EmModel>> &model_ptr);

//    EmAlgorithm() ;

//    EmAlgorithm(int category_count);

    virtual ~EmAlgorithm() {
    }

    virtual std::vector<double> GetParameters();


    virtual void Run() = 0;

    vector<double> GetProportion();

    void PrintSummary();
protected:


    size_t num_category;
    size_t site_count;
    size_t em_count;

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

    vector<vector<double>> temp_stats;


    void Init();


    void ExpectationStep ();

    void MaximizationStep();

    void CalculateProportion();

    virtual void InitialiseProportion();



    virtual void InitialiseParameters() = 0;

    virtual void InitialiseSummaryStat() = 0;

    void ExpectationStep2();


    void ExpectationStep_Old();


    virtual void ExpectationStepCustom(size_t data_index, size_t category_index, double &sum_prob, vector<double> &temp_stat);

    bool EmStoppingCriteria();

    std::vector<double> cache_parameters;

};


static const double EM_CONVERGE_THRESHOLD = 1e-10;
#endif //EM_ALGORITHM_H_
