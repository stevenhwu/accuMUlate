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
    EmSummaryStat(double stat) : stat(stat) {
    }

    EmSummaryStat() {
    }

    virtual ~EmSummaryStat() {
    }

    virtual void print() {
        std::cout << "EMSumStat: " << stat << std::endl;
    }
    virtual void Reset(){
        stat = 0;
    };

    virtual void SetStats(std::vector<double> stats){
        std::cout << "SetStat in EMSumStat: " << stat << std::endl;
        stat = stats[0];
    }
//protprivate:
    const int stat_count = 1;
    double stat = 0;


};


#endif //EM_SUMMARY_STAT_H_
