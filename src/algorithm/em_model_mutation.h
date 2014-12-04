/*
 * em_model_mutation.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_MODEL_MUTATION_H_
#define EM_MODEL_MUTATION_H_


#include "em_model.h"
#include "em_model_mutation.h"

class EmModelMutation : public EmModel {

public:
    EmModelMutation() {
    }

    virtual ~EmModelMutation() {
    }

    virtual void UpdateParameter(double param);



};


#endif //EM_MODEL_MUTATION_H_
