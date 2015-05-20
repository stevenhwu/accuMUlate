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
    cache_log_likelihood = std::numeric_limits<double>::lowest();




    InitialiseProportion();
    InitialiseParameters();
    InitialiseSummaryStat();

    num_blocks = 1;
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


void EmAlgorithmThreadMutation::InitialiseParameters() {
    double low_rate = 1e-10;
    double high_rate = 0.9;



    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> rand_real(1,9);
    std::uniform_int_distribution<> lower(3, 10);

    const int lower_power = lower(gen);
    low_rate = rand_real(gen)*pow(10, -lower_power);

    do {
        std::uniform_int_distribution<> upper(1, lower_power);
        high_rate = rand_real(gen) * pow(10, -upper(gen));
    } while (high_rate < low_rate);

//    low_rate = 1e-10;
//    high_rate = 0.9;
//    low_rate = 1- 1e-10;
//    high_rate = 1- 0.9;
//    model_multi.U

//EM Summary: Ln:-1405.19959316
//Parameters: 3.7007533031e-01	6.5633711554e-09
//Proportions: 1.4027775284e-02	9.8597222472e-01	8.8364417124e-01	5.1913335718e-01	9.8597221824e+01	6.4713016597e-07
//================= END SUMMARY ============


    high_rate =  3.7007533031e-01;
    low_rate = 6.5633711554e-09;
    parameters = {high_rate, low_rate};
    cache_parameters = {high_rate, low_rate};
//    parameters = {low_rate, high_rate};
//    cache_parameters = {low_rate, high_rate};

//FinalSummary: 2.53325e-06	7.12781e-01	9.99326e-01	6.74226e-04	-1471693.98632	7.80611e-11
//FinalSummary: 7.12781e-01	2.53325e-06	6.74226e-04	9.99326e-01	-1471693.98632	7.81574e-11



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
