/*
 * em_algorithm.h
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */


#pragma once

#ifndef EM_ALGORITHM_H_
#define EM_ALGORITHM_H_
#include <iostream>
#include <memory>
#include <vector>
#include <stddef.h>
#include <iostream>
#include <fstream>
#include <atomic>
#include <Eigen/Dense>
#include "em_model.h"
#include "em_data.h"
#include "em_model_mutation.h"
#include "em_logger.h"
#include <mutations/mutation_model_multi_categories.h>

extern const double EM_CONVERGE_THRESHOLD;

/*
Try to be generic, but doesn't work out too well at moment.
Approach 1: completely separate data and model. 1 copy each.
 - repeat data leads to repeated calculation, hard to reuse data
 - hard to pass things around (generically). either non-generic method to pass the parameters or use casting(bad!?)

Approach 2: model incluede data
 - What happen with huge data and huge number of cagetoires. n category = n models

try?? m x d (category x datapoint) strcuture, to avoid any duplications
*/
class EmAlgorithm {


public:

    EmAlgorithm(int num_category0);

    EmAlgorithm(int num_category0, std::vector<std::unique_ptr<EmData>> &data_ptr, EmModel &em_model0);

    EmAlgorithm(std::vector<std::unique_ptr<EmData>> &data_ptr, std::vector<std::unique_ptr<EmModel>> &model_ptr);

    EmAlgorithm(std::vector<std::unique_ptr<EmModel>> &model_ptr);

    virtual ~EmAlgorithm() {
    }

    void Run();

    virtual void RunEM() = 0;

    virtual std::vector<double> GetProportion();

    virtual std::vector<double> GetParameters();

    void PrintSummary();

    void SetOutfilePrefix(const std::string &infile);

    std::string GetEMSummary();

protected:


    size_t num_category;
    size_t site_count;
    size_t max_ite_count;

    size_t stat_count;

    double sum_ratio;
    std::atomic<double> log_likelihood;
//    double log_likelihood;

    EmLogger em_logger;

    std::vector<std::unique_ptr<EmData>> *em_data_ptr;
    std::vector<std::unique_ptr<EmModel>> *em_model_ptr;
//    std::vector<std::shared_ptr<EmModel>> *em_model_ptr;

    std::vector<std::unique_ptr<EmModel>> em_model;
    std::vector<std::unique_ptr<EmSummaryStat>> all_em_stats;

    EmModel *em_model0;//Should be able to remove as well. try to finialse V1 setup

    std::vector<double> parameters;
    std::vector<double> proportion;
    Eigen::ArrayXXd all_probs;

    std::vector<std::vector<double>> temp_stats;
    std::vector<double> cache_parameters;



    void InitWithData();

    void InitWithModel();

    void UpdateEmParameters();

    void ExpectationStepModel();

    void ExpectationStepModelPtr();

    void MaximizationStep();

    void CalculateProportion();

    virtual void InitialiseProportion();

    virtual void InitialiseParameters();

    virtual void InitialiseSummaryStat() = 0;

    virtual void ExpectationStepCustom(size_t data_index, size_t category_index,
            double &sum_prob, std::vector<double> &temp_stat) = 0;

    bool EmStoppingCriteria(int ite);


    void LogEmSummary(int ite);


    std::atomic<double> cache_log_likelihood;

    std::string header;
};


#endif //EM_ALGORITHM_H_
