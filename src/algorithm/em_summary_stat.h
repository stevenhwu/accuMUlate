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
    virtual ~EmSummaryStat() {}


    virtual void Reset() final;
    virtual double GetStat(int index) final;
    virtual size_t GetStatCount() final;

    virtual double MaximiseStats();

    virtual void Print();

    virtual void SetStat(int index, double stat0);
    virtual void SetStats(std::vector<double> stats);

    virtual void UpdateSumWithProportion(double &d, std::unique_ptr<EmSummaryStat> &em_stat_local);
    virtual void UpdateSumWithProportion(double proportion, std::vector<double> &temp_stats);


protected:
    int const stat_count;//int const stat_count;

//    std::vector<double> stat;
//    std::mutex stat_mutex;
//    std::atomic<std::vector<double>> stat;
    std::vector<std::atomic<double>> stat;
};


#endif //EM_SUMMARY_STAT_H_
