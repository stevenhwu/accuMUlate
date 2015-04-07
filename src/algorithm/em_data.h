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
#include "em_summary_stat.h"
#include "em_model.h"
#include "em_model_mutation_v1.h"
#include "em_model_binomial.h"
/*
EmData was used with original version of EmModel, where model doesn't contains data,
i.e. keep model and data completely separate
This turns out might not be a good idea
1) many identical data needs to be recalculated
2) hard to be generic, might (at least to my knowledge) needs some sort of casting
3) some strange pointers usage.

maybe removed the whole structure later.
It's is currently used for testing the original implementations
*/
class EmData {


public:

    virtual ~EmData() {
    }
    virtual void UpdateSummaryStat(double &prob, std::vector<double> &temp_stat) = 0;
    virtual void UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat) = 0;
//    virtual void UpdateEmModel(EmModel &em_model) = 0;

    virtual void UpdateEmModel(EmModel *em_model) = 0;

//    virtual std::unordered_map<std::string ,std::vector<double>> GetDataVector() const;
//    virtual void UpdateEmModel(std::unique_ptr<EmModel> &em_model)=0;
//    virtual void UpdateEmModel(std::unique_ptr<EmModelMutation> &em_model) = 0;
};


#endif //EM_DADA_H_
