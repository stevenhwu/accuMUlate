/*
 * em_algorithm_multithread.cc
 *
 *  Created on: 4/21/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/21/15.
//

#include <atomic>
#include "em_algorithm_thread.h"

EmAlgorithmMultiThreading::EmAlgorithmMultiThreading(MutationModelMultiCategories &model_multi0, uint32_t thread_count0)
        : EmAlgorithm(model_multi0.GetCategoriesCount()), model_multi(model_multi0), thread_count(thread_count0) {
    std::cout << "ModelMultiCategories MultiThread Construct: " << num_category << "\t" << "Thread_count: " << "\t" <<
    thread_count << "\t One data, one model" << std::endl;
}


void EmAlgorithmMultiThreading::ExpectationStepModelPtrMTMulti() {
    //Non generic method fro em_algorithm_thread_mutation ONLY

//    UpdateEmParameters();
    for (size_t r = 0; r < num_category; ++r) {
        all_em_stats[r]->Reset();
        model_multi.UpdateOneMinusExpBeta(r, parameters[r]); //1-exp_beta
    }
    log_likelihood = 0;

    {
        ThreadPool tp(thread_count);
//        for (int i = 0; i < num_blocks - 1; ++i) {
//            tp.enqueue(boost::bind(&EmAlgorithmMultiThreading::MultiCategories, this, i * block_size, (i + 1) * block_size));
//        }
//        tp.enqueue(boost::bind(&EmAlgorithmMultiThreading::MultiCategories, this, (num_blocks - 1) * block_size, site_count));
        for (int i = 0; i < num_blocks ; ++i) {
            tp.enqueue(boost::bind(&EmAlgorithmMultiThreading::MultiCategories, this, blocks_info[i].first, blocks_info[i].second));
        }

    }


}


void EmAlgorithmMultiThreading::MultiCategories(size_t site_start, size_t site_end) {

//    std::cout << "GetID: " << boost::this_thread::get_id() << std::dec<<  "\t" << ((int) site_start )<< "\t" << (int) site_end  << std::endl;

    double thread_likelihood = 0;
    std::vector<std::vector<double>> thread_stats(num_category, std::vector<double>(stat_count, 0));

    double log_likelihood_scaler = 0;
    double sum_prob = 0;
    double stat_diff = 0;
//    std::vector<std::vector<double>> temp_stats (num_category, std::vector<double>(stat_count, 0));
    std::vector<double> temp_stats (num_category, 0);
    std::vector<double> local_probs (num_category, 0);

    for (size_t site = site_start; site < site_end; ++site) {
        int descendant_count = model_multi.GetDescendantCount(site);
        double sum = 0;
        for (size_t r = 0; r < num_category; ++r) {
            model_multi.CalculateAncestorToDescendant(r, site, sum_prob, stat_diff, log_likelihood_scaler);
            temp_stats[r] = stat_diff;

            local_probs[r] = proportion[r] * sum_prob;
            sum += local_probs[r];
        }

        thread_likelihood += log(sum) + log_likelihood_scaler;

        for (size_t r = 0; r < num_category; ++r) {
            double prob = local_probs[r] / sum;

            thread_stats[r][1] += prob * temp_stats[r];
            thread_stats[r][0] += prob * (descendant_count-temp_stats[r]);
//            std::cout << descendant_count << "\t" << temp_stats[r] << "\t" << prob << "\t" << (prob*temp_stats[r])<< std::endl;
            all_probs(r, site) = local_probs[r];//threads
        }

    }

    auto current = log_likelihood.load();
    while (!log_likelihood.compare_exchange_weak(current, current + thread_likelihood));//threads
    for (size_t r = 0; r < num_category; ++r) {
        all_em_stats[r]->UpdateSumWithProportionSynchronized(thread_stats[r]);//threads
    }

}


////////////////////////////////////////////////////////////////////////////////

void EmAlgorithmMultiThreading::ExpectationStepModelPtrMT() {
//For em_algorithm_mutation.cc Old code

//    std::cout << "GetID: " << boost::this_thread::get_id() << '\n';
//    std::cout << boost::thread::hardware_concurrency() << '\n';
//
////    boost::scoped_thread<> t{boost::thread{thread}};
//    boost::thread t1{fill}, t2{print}, t3{count};
//    t1.join();
//    t2.join();
//    t3.join();
//    std::cout << "Sum: " << sum << '\n';

//    size_t block_size = 5000;

//    const int thread_count = 4;

    char temp[1000];
    UpdateEmParameters();
//
//    boost::thread t[thread_count];
    log_likelihood = 0;
//    std::cout << thread_count << "\t" << block_size<< std::endl;

//    size_t block_size = site_count/thread_count;
//    for (int i = 0; i < thread_count-1; ++i) {
//        t[i] = boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, i * block_size, (i + 1) * block_size);
//    }
//    t[thread_count-1] = boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, (thread_count-1)*block_size,site_count);
////    std::cout << "==All Ln: " << log_likelihood << std::endl;
//    for (int i = 0; i < thread_count; ++i)
//        t[i].join();

//    size_t block_size = site_count/thread_count;
//    thread_vector.clear();
//    for (int i = 0; i < thread_count-1; ++i) {
//        thread_vector.push_back(
//                boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, i * block_size, (i + 1) * block_size));
//    }
//    thread_vector.push_back(boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, (thread_count-1)*block_size,site_count));
//
//    for (int i = 0; i < thread_count; ++i){
//        thread_vector[i].join();
//    }



//    boost::asio::io_service io_service;
//    boost::asio::io_service::work work(io_service);
//    boost::thread_group threads;
//    int thread_count  = 2;
//    for (std::size_t i = 0; i < thread_count; ++i) {
//        threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
//    }
    size_t num_blocks = 1000;
    size_t block_size = site_count / num_blocks;
//     for (int i = 0; i < num_blocks-1; ++i) {
//         std::cout << i << "\t" << block_size << std::endl;
//         io_service.post(boost::bind(&EmAlgorithmMultiThreading::WorkingThread, this, i * block_size, (i + 1) * block_size)  );
//     }
//    io_service.post(boost::bind(&EmAlgorithmMultiThreading::WorkingThread, this, (num_blocks-1)*block_size,site_count) );
//    io_service.stop();
//    threads.join_all();

    {
        ThreadPool tp(1);
        for (int i = 0; i < num_blocks - 1; ++i) {
            tp.enqueue(boost::bind(&EmAlgorithmMultiThreading::WorkingThread, this, i * block_size, (i + 1) * block_size));
        }
        tp.enqueue(boost::bind(&EmAlgorithmMultiThreading::WorkingThread, this, (num_blocks - 1) * block_size, site_count));

//        for (int i = 0; i < num_blocks - 1; ++i) {
//            tp.enqueue(boost::bind(&EmAlgorithmMultiThreading::EmptyThread, this, i * block_size, (i + 1) * block_size));
//        }
//        tp.enqueue(boost::bind(&EmAlgorithmMultiThreading::EmptyThread, this, (num_blocks - 1) * block_size, site_count));
    }


//    boost::this_thread::sleep_for(boost::chrono::seconds{3});


//    sprintf(temp, "%.15f\t", log_likelihood.load());
////    sprintf(temp, "%.15f\t", log_likelihood);
//    std::cout << "\n==MT==All Ln: " << temp << std::endl;
//    for (size_t r = 0; r < num_category; ++r) {
//        all_em_stats[r]->Print();
//    }





//    for (size_t r = 0; r < num_category; ++r) {
//        all_em_stats[r]->Print();
//    }

//    UpdateEmParameters();
//    log_likelihood = 0;
//    WorkingThread(0,site_count);
//////    WorkingThread(0,31321);
//
//    sprintf(temp, "%.15f\t", log_likelihood.load());
//    std::cout << "\n==Single==All Ln: " << temp << std::endl;
//    for (size_t r = 0; r < num_category; ++r) {
//        all_em_stats[r]->Print();
//    }
//    std::exit(12);
//7.204842241669e-01	9.095599810937e-07	1.305992121779e-04	9.998694007878e-01	9.998694007878e-01	Time new: 61	61484826

}
//__attribute_deprecated__
[[deprecated]]
void EmAlgorithmMultiThreading::WorkingThread(size_t site_start, size_t site_end) {

    double log_likelihood_scaler = 0;
    double temp_likelihood = 0;
    int stat_count = 2;
    std::vector<std::vector<double>> temp_stats (num_category, std::vector<double>(2, 0));
    std::vector<EmSummaryStat> local_stat (num_category, stat_count);


    for (size_t s = site_start; s < site_end; ++s) {

        double sum_prob = 0;
        for (size_t r = 0; r < num_category; ++r) {
            (*em_model_ptr)[r]->UpdateSummaryStat(s, sum_prob, temp_stats[r], log_likelihood_scaler);
            all_probs(r, s) = proportion[r] * sum_prob;
        }
        temp_likelihood += log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler;

        double sum = all_probs.col(s).sum();
        for (size_t r = 0; r < num_category; ++r) {
            double prob = all_probs(r,s) / sum;

//            all_em_stats[r]->UpdateSumWithProportion(prob, temp_stats[r]);
            local_stat[r].UpdateSumWithProportion(prob, temp_stats[r]);
        }
    }
    auto current = log_likelihood.load();
    while (!log_likelihood.compare_exchange_weak(current, current + temp_likelihood));
    for (size_t r = 0; r < num_category; ++r) {
        all_em_stats[r]->UpdateSumWithProportionSynchronized(local_stat[r].GetStats());
    }

//    for (size_t r = 0; r < num_category; ++r) {
//        all_em_stats[r]->Print();
//    }
//    std::cout << "====Ln: " << log_likelihood << std::endl;

////    printf("Ln: %.40f\n", log_likelihood);

}
