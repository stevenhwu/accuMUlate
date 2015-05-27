#include "em_algorithm_thread_mutation.h"
#include "em_summary_stat_mutation.h"

EmAlgorithmThreadMutation::EmAlgorithmThreadMutation(MutationModelMultiCategories &model_multi0, uint32_t thread_count0)
        : EmAlgorithmMultiThreading(model_multi0, thread_count0) {

    site_count = model_multi.GetSiteCount();
    stat_count = EmSummaryStatMutation::EM_SUMMARY_STAT_MUTATION_STATS_COUNT;
    std::cout << "EM Site_count: " << site_count << "\t" << "Thread_count" << "\t" << thread_count << std::endl;
    all_probs = Eigen::ArrayXXd::Zero(num_category, site_count);
    parameters = std::vector<double>(num_category, 0);
    cache_parameters = std::vector<double>(num_category, 0);
    log_likelihood = 0;//std::numeric_limits<double>::lowest();
    cache_log_likelihood = std::numeric_limits<double>::lowest();




    InitialiseProportion();
    InitialiseParameters();
    InitialiseSummaryStat();

    num_blocks = 120;
    blocks_info.reserve(num_blocks);
    size_t block_size = site_count / num_blocks;
    for (int i = 0; i < num_blocks - 1; ++i) {
        blocks_info.emplace_back(i * block_size, (i + 1) * block_size);
    }
    blocks_info.emplace_back( (num_blocks - 1) * block_size, site_count);


}

EmAlgorithmThreadMutation::~EmAlgorithmThreadMutation() {

}



void EmAlgorithmThreadMutation::RunEM() {

    size_t i = 0;
    bool notConverged = true;
    while(notConverged){
        ExpectationStepModelPtrMTMulti();
        MaximizationStep();
        notConverged = EmStoppingCriteria(i);
        i++;

    }

}




void EmAlgorithmThreadMutation::InitialiseSummaryStat() {

//    stat_count = EmSummaryStatMutation::EM_SUMMARY_STAT_MUTATION_STATS_COUNT;
//    temp_stats = std::vector<std::vector<double>>(num_category);
    for (size_t i = 0; i < num_category; ++i) {
        all_em_stats.emplace_back(new EmSummaryStatMutation());
//        temp_stats[i] = std::vector<double>(all_em_stats[i]->GetStatCount());
    }


}


void EmAlgorithmThreadMutation::ExpectationStepCustom(size_t data_index, size_t category_index,
        double &sum_prob, std::vector<double> &temp_stat) {
    std::cout << "Error!! should NOT call ExpectationStepCustom here" << std::endl;
    exit(40);
}
