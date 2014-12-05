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


public:

    EmSummaryStat();
    virtual ~EmSummaryStat() {
    }


    double GetStat(int index);

    void SetStat(int index, double stat0);


    virtual void print();
    virtual void Reset();

    virtual void UpdateSumWithProportion(double &d, EmSummaryStat &mutation);

    virtual void SetStats(std::vector<double> stats);

protected:
    int stat_count = 1;
    std::vector<double> stat;


};


#endif //EM_SUMMARY_STAT_H_
