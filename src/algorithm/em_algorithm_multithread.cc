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

    const int num_thread = 4;
    size_t block_size = site_count/num_thread;

    UpdateEmParameters();
//
//    boost::thread t[num_thread];
    log_likelihood = 0;
//    std::cout << num_thread << "\t" << block_size<< std::endl;
    for (int i = 0; i < num_thread-1; ++i) {
        t[i] = boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, i * block_size, (i + 1) * block_size);
    }
    t[num_thread-1] = boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, (num_thread-1)*block_size,site_count);
//    t[0] = boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, 0,31321);
//    t[0] = boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, 0,15660);
//    t[1] = boost::thread(&EmAlgorithmMultiThreading::WorkingThread, this, 15660, 31321);
//
//    std::cout << "==All Ln: " << log_likelihood << std::endl;
    for (int i = 0; i < num_thread; ++i)
        t[i].join();


//    for (size_t r = 0; r < num_category; ++r) {
//        all_em_stats[r]->Print();
//    }

//    UpdateEmParameters();
//    log_likelihood = 0;
//    WorkingThread(0,num_thread*block_size);
////    WorkingThread(0,31321);
//    std::cout << "\n\nSingle\n==All Ln: " << log_likelihood << std::endl;
//    for (size_t r = 0; r < num_category; ++r) {
//        all_em_stats[r]->Print();
//    }
//    std::exit(12);

}

//static void WorkingThread(size_t site_start, size_t site_end);

void EmAlgorithmMultiThreading::WorkingThread(size_t site_start, size_t site_end) {

//    std::cout << "GetID: " << boost::this_thread::get_id() << "\t" << site_start << "\t" << site_end << '\n';

    double log_likelihood_scaler = 0;
//    std::vector<std::vector<double>> temp_stats (2);
//    for (size_t i = 0; i < num_category; ++i) {
//        std::cout << i << std::endl;
//        temp_stats[i] = std::vector<double>(2);
//        temp_stats[i].assign(2,0);
//        std::cout << temp_stats[i][0] << "\t" << temp_stats[i][1] << std::endl;
//    }
//    for (size_t s = 0; s < site_count; ++s) {
    for (size_t s = site_start; s < site_end; ++s) {
        double sum_prob = 0;
        for (size_t r = 0; r < num_category; ++r) {
            (*em_model_ptr)[r]->UpdateSummaryStat(s, sum_prob, temp_stats[r], log_likelihood_scaler);
            all_probs(r, s) = proportion[r] * sum_prob;
        }

//        log_likelihood = log_likelihood + log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler  ;
        double temp = log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler;
        auto current = log_likelihood.load();
        while (!log_likelihood.compare_exchange_weak(current, current + temp));

        double sum = all_probs.col(s).sum();
//        auto a = all_probs.col(s);
//        auto aa = a.sum();

//        std::cout << "P: "<<all_probs(0,s) << "\t" << all_probs(1,s)<<  "\t" << sum << std::endl;
//        std::cout << a(0) << "\t" << a[1] << "\t" << aa << std::endl;
        for (size_t r = 0; r < num_category; ++r) {
            double prob = all_probs(r,s) / sum;
//            std::cout << all_probs(r,s) << "\t" <<  sum << std::endl;
            all_em_stats[r]->UpdateSumWithProportion(prob, temp_stats[r]);
        }
    }
//    for (size_t r = 0; r < num_category; ++r) {
//        all_em_stats[r]->Print();
//    }
//    std::cout << "====Ln: " << log_likelihood << std::endl;
////    printf("Ln: %.40f\n", log_likelihood);

}