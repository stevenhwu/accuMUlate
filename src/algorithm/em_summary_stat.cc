#include "em_summary_stat.h"

EmSummaryStat::EmSummaryStat() : EmSummaryStat(EM_SUMMARY_STAT_STATS_COUNT){
};


EmSummaryStat::EmSummaryStat(int const stat_count0) : stat_count(stat_count0) {
//    stat.assign(stat_count, std::atomic<double> {0});
//    stat.resize(stat_count);
//    std::cout << "creat stat:" << "\t" << stat_count << std::endl;
    stat.assign(stat_count, 0);
//    for (int i = 0; i < stat_count; ++i) {
//        stat.emplace_back(0);
//    }
}

size_t EmSummaryStat::GetStatCount(){
    return stat_count;
}

void EmSummaryStat::Print() {
    std::cout << "Print_EmSumStat: C=" << stat_count <<" ";
    for (int i = 0; i < stat_count; ++i) {
        std::cout << stat[i] << "\t" ;
    }
    std::cout << std::endl;
}


double EmSummaryStat::GetStat(int index) {
    return stat[index];
}

std::vector<double>& EmSummaryStat::GetStats() {
    return stat;
}

void EmSummaryStat::Reset() {
    for (int i = 0; i < stat_count; ++i) {
        stat[i] = 0;
    }
}



void EmSummaryStat::UpdateSumWithProportion(double &d, std::unique_ptr<EmSummaryStat> &summary_stat) {
//    for (int i = 0; i < summary_stat->stat_count; ++i) {
//        stat[i] += d * summary_stat->GetStat(i);
//        double temp = d * summary_stat->GetStat(i);
//        std::cout << "VECTOR STAT: " << stat.size() << "\t" << i<<"\t" << temp << "\t" << stat[i] <<std::endl;
//    }
    for (int i = 0; i < summary_stat->stat_count; ++i) {
//        stat[i] += d * summary_stat->GetStat(i);
        double temp = d * summary_stat->GetStat(i);
        std::cout << "VECTOR STAT: " << stat.size() << "\t" << i<<"\t" << d << "\t" <<  summary_stat->GetStat(i) << "\t" <<


            temp << "\t" << stat[i] <<std::endl;
    }
}

void EmSummaryStat::SetStat(int index, double stat0) {
    stat[index] = stat0;
}

void EmSummaryStat::SetStats(std::vector<double> stats0) {
    for (int i = 0; i < stat_count; ++i) {
        stat[i] = stats0[i];
    }

}

double EmSummaryStat::MaximiseStats() {
    //TODO: What is the default here?? Or make this a abstract class
    exit(-100);
}

void EmSummaryStat::UpdateSumWithProportion(double proportion, std::vector<double> &temp_stats) {
//    stat_mutex.lock();
    for (int i = 0; i < stat_count; ++i) {
//        std::cout << proportion << "\t" << temp_stats[i] << std::endl;
        stat[i] += proportion * temp_stats[i];

//        auto current = stat[i].load();
//        auto temp = proportion * temp_stats[i];
//        while (!stat[i].compare_exchange_weak(current, current + temp));

    }
//    stat_mutex.unlock();
}

void EmSummaryStat::UpdateSumWithProportionSynchronized(std::vector<double> &temp_stats) {
    std::lock_guard<std::mutex> lock(stat_mutex);
    for (int i = 0; i < stat_count; ++i) {
        stat[i] += temp_stats[i];
    }
}

void EmSummaryStat::ChangeStatCount(int new_count) {
    stat_count = new_count;
    stat.assign(stat_count, 0);

}
