/*
 * em_model.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_MODEL_H_
#define EM_MODEL_H_


class EmModel {
public:
    virtual void UpdateParameter(double param) = 0;
};


#endif //EM_MODEL_H_
