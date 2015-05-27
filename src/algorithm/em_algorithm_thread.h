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

#include <boost/asio/io_service.hpp>


#include "em_algorithm.h"


typedef std::unique_ptr<boost::asio::io_service::work> asio_worker;

struct ThreadPool {
    ThreadPool(size_t threads) :service(), working(new asio_worker::element_type(service)) {
        while(threads--)
        {
            auto worker = boost::bind(&boost::asio::io_service::run, &(this->service));
            g.add_thread(new boost::thread(worker));
        }
    }

    template<class F>
    void enqueue(F f){
        service.post(f);
    }

    ~ThreadPool() {
        working.reset(); //allow run() to exit
        g.join_all();
        service.stop();
    }

private:
    boost::asio::io_service service; //< the io_service we are wrapping
    asio_worker working;
    boost::thread_group g; //< need to keep track of threads so we can join them
};


class EmAlgorithmMultiThreading : public EmAlgorithm {

public:

    EmAlgorithmMultiThreading(MutationModelMultiCategories &model_multi0, uint32_t thread_count0);


    void ExpectationStepModelPtrMTMulti();

    void MultiCategories(size_t site_start, size_t site_end);



    void ExpectationStepModelPtrMT();

    [[deprecated]]
    void WorkingThread(size_t site_start, size_t site_end);



protected:

    MutationModelMultiCategories model_multi;
    std::vector<std::pair<size_t, size_t>> blocks_info;


    uint32_t thread_count;
    size_t num_blocks;

    std::mutex lock_stat;
};


#endif //EM_ALGORITHM_MULTITHREAD_H

