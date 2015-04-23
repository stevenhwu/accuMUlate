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
#include "em_algorithm_multithread.h"




void thread()
{
    for (int i = 0; i < 5; ++i)
    {
//        wait(0.1);
        std::cout << i << '\n';
    }
}
void wait(int seconds)
{
    boost::this_thread::sleep_for(boost::chrono::seconds{seconds});
}

boost::shared_mutex mutex;
std::vector<int> random_numbers;

void fill()
{
    std::srand(static_cast<unsigned int>(std::time(0)));
    for (int i = 0; i < 3; ++i)
    {
        boost::unique_lock<boost::shared_mutex> lock{mutex};
        random_numbers.push_back(std::rand());
        lock.unlock();
        wait(1);
    }
}

void print()
{
    for (int i = 0; i < 3; ++i)
    {
        wait(1);
        boost::shared_lock<boost::shared_mutex> lock{mutex};
        std::cout << random_numbers.back() << '\n';
    }
}

int sum = 0;
boost::mutex mutex2;
void count()
{
    for (int i = 0; i < 3; ++i)
    {
        wait(1);
        boost::shared_lock<boost::shared_mutex> lock{mutex};
        sum += random_numbers.back();
    }
}

void init()
{
    static boost::thread_specific_ptr<bool> tls;
    if (!tls.get())
    {
        tls.reset(new bool{true});
        boost::lock_guard<boost::mutex> lock{mutex2};
        std::cout << "done" << '\n';
    }
}

void thread2()
{
    init();
    init();
}



void EmAlgorithmMultiThreading::ExpectationStepModelPtrMT() {


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

//    const int num_thread = 4;

    char temp[1000];
    UpdateEmParameters();
//
//    boost::thread t[num_thread];
    log_likelihood = 0;
//    std::cout << num_thread << "\t" << block_size<< std::endl;

//    size_t block_size = site_count/num_thread;
//    for (int i = 0; i < num_thread-1; ++i) {
//        t[i] = boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, i * block_size, (i + 1) * block_size);
//    }
//    t[num_thread-1] = boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, (num_thread-1)*block_size,site_count);
////    std::cout << "==All Ln: " << log_likelihood << std::endl;
//    for (int i = 0; i < num_thread; ++i)
//        t[i].join();

//    size_t block_size = site_count/num_thread;
//    thread_vector.clear();
//    for (int i = 0; i < num_thread-1; ++i) {
//        thread_vector.push_back(
//                boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, i * block_size, (i + 1) * block_size));
//    }
//    thread_vector.push_back(boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, (num_thread-1)*block_size,site_count));
//
//    for (int i = 0; i < num_thread; ++i){
//        thread_vector[i].join();
//    }



//    boost::asio::io_service io_service;
//    boost::asio::io_service::work work(io_service);
//    boost::thread_group threads;
//    int num_thread  = 2;
//    for (std::size_t i = 0; i < num_thread; ++i) {
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
void EmAlgorithmMultiThreading::EmptyThread(size_t site_start, size_t site_end) {
//    double temp = site_start + site_end;
//    auto current = log_likelihood.load();
//    while (!log_likelihood.compare_exchange_weak(current, current + temp));
//
//    std::vector<double> temp_v;
//    temp_v.push_back(site_start);
//    temp_v.push_back(site_end);
//    for (size_t r = 0; r < num_category; ++r) {
////        double prob = all_probs(r,s) / sum;
////            std::cout << all_probs(r,s) << "\t" <<  sum << std::endl;
//
//        all_em_stats[r]->UpdateSumWithProportion( (double) (site_start/(site_end*1.0)), temp_v);
//    }
//    for (size_t s = site_start; s < site_end; ++s) {
////        std::cout << s << std::endl;
//
//        for (size_t r = 0; r < num_category; ++r) {
//
//            all_probs(r, s) = (site_start/(site_end*1.0))*(1.0/s);
//        }
//    }

    double log_likelihood_scaler = 0;
    std::vector<std::vector<double>> temp_stats (2);
    for (size_t i = 0; i < num_category; ++i) {
        temp_stats[i] = std::vector<double>(2);
//        temp_stats[i].assign(2,0);
    }

//    for (size_t s = 0; s < site_count; ++s) {
    double temp_likelihood = 0;
    for (size_t s = site_start; s < site_end; ++s) {
//        std::cout << s << std::endl;
        double sum_prob = 1.0;//*s*site_start/site_end;;
        for (size_t r = 0; r < 2; ++r) {
            (*em_model_ptr)[r]->UpdateSummaryStat(s, sum_prob, temp_stats[r], log_likelihood_scaler);
            all_probs(r, s) = proportion[r] * sum_prob;
        }
//        pModel0->UpdateSummaryStat(s, sum_prob, temp_stats[0], log_likelihood_scaler);
//        lock.lock();// 142 vs 40
//        pModel1->UpdateSummaryStat(s, sum_prob, temp_stats[1], log_likelihood_scaler);
//    lock.unlock();
//        log_likelihood = log_likelihood + log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler  ;
        temp_likelihood += log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler;
//        double temp = log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler  ;
//        auto current = log_likelihood.load();
//        while (!log_likelihood.compare_exchange_weak(current, current + temp));

        double sum = all_probs.col(s).sum();

        for (size_t r = 0; r < num_category; ++r) {
            double prob = all_probs(r,s) / sum;

            all_em_stats[r]->UpdateSumWithProportion(prob, temp_stats[r]);
        }
    }

    auto current = log_likelihood.load();
    while (!log_likelihood.compare_exchange_weak(current, current + temp_likelihood));

}
void EmAlgorithmMultiThreading::WorkingThread(size_t site_start, size_t site_end) {

//    std::cout << "GetID: " << boost::this_thread::get_id() << std::dec<<  "\t" << ((int) site_start )<< "\t" << (int) site_end  << std::endl;

    double log_likelihood_scaler = 0;
    double temp_likelihood = 0;
    int stat_count = 2;
    std::vector<std::vector<double>> temp_stats (num_category, std::vector<double>(2, 0));
    std::vector<EmSummaryStat> local_stat (num_category, stat_count);



//    for (size_t s = 0; s < site_count; ++s) {
    for (size_t s = site_start; s < site_end; ++s) {
//        std::cout << s << std::endl;
        double sum_prob = 0;
        for (size_t r = 0; r < num_category; ++r) {
            (*em_model_ptr)[r]->UpdateSummaryStat(s, sum_prob, temp_stats[r], log_likelihood_scaler);
            all_probs(r, s) = proportion[r] * sum_prob;
        }
//        std::exit(3);

//        log_likelihood = log_likelihood + log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler  ;
//        double temp = log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler;
//        auto current = log_likelihood.load();
//        while (!log_likelihood.compare_exchange_weak(current, current + temp));
        temp_likelihood += log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler;

        double sum = all_probs.col(s).sum();
//        auto a = all_probs.col(s);
//        auto aa = a.sum();
//        std::cout << "P: "<<all_probs(0,s) << "\t" << all_probs(1,s)<<  "\t" << sum << std::endl;
//        std::cout << a(0) << "\t" << a[1] << "\t" << aa << std::endl;
        for (size_t r = 0; r < num_category; ++r) {
            double prob = all_probs(r,s) / sum;
//            std::cout << all_probs(r,s) << "\t" <<  sum << std::endl;
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

void EmAlgorithmMultiThreading::ExpectationStepModelPtrMTMulti() {
    char temp[1000];
//    UpdateEmParameters();
    for (size_t r = 0; r < num_category; ++r) {
        all_em_stats[r]->Reset();
        model_multi.UpdateExpBeta(r,parameters[r]); //exp_beta
    }
    log_likelihood = 0;
    size_t num_blocks = 120;
    size_t block_size = site_count / num_blocks;
    {
        ThreadPool tp(3);
        for (int i = 0; i < num_blocks - 1; ++i) {
            tp.enqueue(boost::bind(&EmAlgorithmMultiThreading::MultiCategories, this, i * block_size, (i + 1) * block_size));
        }
        tp.enqueue(boost::bind(&EmAlgorithmMultiThreading::MultiCategories, this, (num_blocks - 1) * block_size, site_count));
    }
//    sprintf(temp, "%.15f\t", log_likelihood.load());
//    sprintf(temp, "%.15f\t", log_likelihood);
//    std::cout << "\n==MultiCategories==All Ln: " << temp << std::endl;
//7.204842241669e-01	9.095599810937e-07	1.305992121779e-04	9.998694007878e-01	9.998694007878e-01	Time new: 61	61484826
}

void EmAlgorithmMultiThreading::MultiCategories(size_t site_start, size_t site_end) {

//    std::cout << "GetID: " << boost::this_thread::get_id() << std::dec<<  "\t" << ((int) site_start )<< "\t" << (int) site_end  << std::endl;

    double log_likelihood_scaler = 0;
    double temp_likelihood = 0;
    int stat_count = 2;

    double prob = 0;
    double stat_diff = 0;

    std::vector<std::vector<double>> temp_stats (num_category, std::vector<double>(2, 0));
    std::vector<std::vector<double>> local_stats (num_category, std::vector<double>(2, 0));
    std::vector<double> local_probs (num_category, 0);

    for (size_t s = site_start; s < site_end; ++s) {
        double sum_prob = 0;
        for (size_t r = 0; r < num_category; ++r) {

            model_multi.CalculateAncestorToDescendant(r, s, prob, stat_diff, log_likelihood_scaler);
            temp_stats[r][0] = 1-stat_diff;
            temp_stats[r][1] = stat_diff;

            local_probs[r] = proportion[r] * prob;
        }
//        std::exit(3);

        double sum = local_probs[0]+local_probs[1];
        temp_likelihood += log(sum) + log_likelihood_scaler;

        for (size_t r = 0; r < num_category; ++r) {
            double prob = local_probs[r] / sum;
            local_stats[r][0] += prob * temp_stats[r][0];
            local_stats[r][1] += prob * temp_stats[r][1];
        }
        all_probs(0,s) = local_probs[0];
        all_probs(1,s) = local_probs[1];

    }
    auto current = log_likelihood.load();
    while (!log_likelihood.compare_exchange_weak(current, current + temp_likelihood));
    for (size_t r = 0; r < num_category; ++r) {
        all_em_stats[r]->UpdateSumWithProportionSynchronized(local_stats[r]);
    }

}
