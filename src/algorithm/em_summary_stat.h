/*
 * em_summary_stat.h
 *
 *  Created on: 12/4/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef EM_SUMMARY_STAT_H_
#define EM_SUMMARY_STAT_H_

#include <vector>
class EmSummaryStat {

    static const int EM_SUMMARY_STAT_STATS_COUNT = 1;
public:
    EmSummaryStat();
    EmSummaryStat(int const stat_count0);
    virtual ~EmSummaryStat() {
    }


    double GetStat(int index);
    //CHECK: These don't have to be virtual right??
    virtual void SetStat(int index, double stat0);
    virtual void SetStats(std::vector<double> stats);

    virtual void print();
    virtual void Reset();
    virtual void UpdateSumWithProportion(double &d, EmSummaryStat &mutation);



    virtual double MaximiseStats();

protected:
    int const stat_count;//int const stat_count;
    std::vector<double> stat;


};


#endif //EM_SUMMARY_STAT_H_
