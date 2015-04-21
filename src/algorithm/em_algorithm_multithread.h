/*
 * em_algorithm_multithread.cc
 *
 *  Created on: 4/21/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/21/15.
//

#pragma once
#ifndef EM_ALGORITHM_MULTITHREAD_H
#define EM_ALGORITHM_MULTITHREAD_H

#include <boost/thread.hpp>
#include <boost/thread/scoped_thread.hpp>
#include <boost/chrono.hpp>

#include "em_algorithm.h"


class EmAlgorithmMultiThreading : public EmAlgorithm {

public:


    EmAlgorithmMultiThreading(std::vector<std::unique_ptr<EmModel>> &model_ptr) : EmAlgorithm(model_ptr) {


    }
    static const int num_thread = 4;
    boost::thread t[num_thread];

    void ExpectationStepModelPtrMT();
    void WorkingThread(size_t site_start, size_t site_end);
};


#endif //EM_ALGORITHM_MULTITHREAD_H

