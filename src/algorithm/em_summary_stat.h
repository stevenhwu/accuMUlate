/*
 * em_summary_stat.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_SUMMARY_STAT_H_
#define EM_SUMMARY_STAT_H_

#include <iostream>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <mutex>
#include <atomic>

//#include <stddef.h>

class EmSummaryStat {

    static const int EM_SUMMARY_STAT_STATS_COUNT = 1;
public:
    EmSummaryStat();
    EmSummaryStat(int const stat_count0);



    EmSummaryStat(const EmSummaryStat& other){
//        this = EmSummaryStat(2);
        this->ChangeStatCount(other.stat_count);
//        std::cout << other.stat_count << "\t" << this->GetStatCount() << "\t" << other.stat_count << "\t" << this->stat_count<< std::endl;
//        std::cout << other.stat[0] << std::endl;
//        std::cout << other.stat[1] << std::endl;
//        for (int i = 0; i < other.stat_count; ++i) {
////            std::cout << "set:" << i << std::endl;
//            this->SetStat(i, other.stat[i]);
//        }
        this->stat = other.stat;
    }

    virtual ~EmSummaryStat() = default;

    virtual void Reset() final;
    virtual double GetStat(int index) final;
    virtual size_t GetStatCount() final;


    virtual std::vector<double>& GetStats() final;
    virtual void SetStat(int index, double stat0) final;
    virtual void SetStats(std::vector<double> stats) final;



    virtual double MaximiseStats();
    virtual void Print();
    virtual void UpdateSumWithProportion(double &d, std::unique_ptr<EmSummaryStat> &em_stat_local);
    virtual void UpdateSumWithProportion(double proportion, std::vector<double> &temp_stats);

    virtual void UpdateSumWithProportionSynchronized(std::vector<double> &temp_stats);


protected:
    int stat_count;

    std::vector<double> stat;
    std::mutex stat_mutex;
//    std::lock_guard<std::mutex> lock(g_i_mutex);

//    std::atomic<std::vector<double>> stat;
//    std::vector<std::atomic<double>> stat ;
//    boost::lockfree::queue<double> stas;

private:
    void ChangeStatCount(int new_count);

};


#endif //EM_SUMMARY_STAT_H_
