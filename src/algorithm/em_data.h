/*
 * em_data.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef EM_DADA_H_
#define EM_DADA_H_


#include <iostream>

//#include <unordered_map>
#include "em_summary_stat.h"
#include "em_model.h"
#include "em_model_mutation_v1.h"
#include "em_model_binomial.h"

class EmData {

public:

    virtual ~EmData() {
    };
    virtual void UpdateSummaryStat(double &prob, std::vector<double> &temp_stat) = 0;
    virtual void UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat) = 0;
//    virtual void UpdateEmModel(EmModel &em_model) = 0;

    virtual void UpdateEmModel(EmModel *em_model) = 0;

//    virtual std::unordered_map<std::string ,std::vector<double>> GetDataVector() const;
//    virtual void UpdateEmModel(std::unique_ptr<EmModel> &em_model)=0;
//    virtual void UpdateEmModel(std::unique_ptr<EmModelMutation> &em_model) = 0;
};


#endif //EM_DADA_H_
