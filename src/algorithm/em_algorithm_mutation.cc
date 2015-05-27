#include "em_algorithm_mutation.h"
#include "em_summary_stat_mutation.h"


EmAlgorithmMutation::EmAlgorithmMutation(std::vector<std::unique_ptr<EmModel>> &model_ptr) : EmAlgorithm(model_ptr) {

    InitWithModel();
    stat_count = 2;
    cache_log_likelihood = std::numeric_limits<double>::lowest();

}


EmAlgorithmMutation::~EmAlgorithmMutation() {

}

void EmAlgorithmMutation::RunEM() {

    size_t i = 0;
    bool notConverged = true;
    while(notConverged){
        ExpectationStepModelPtr();

        MaximizationStep();
        notConverged = EmStoppingCriteria(i);
        i++;
    }

}
//
//void EmAlgorithmMutation::InitialiseParameters() {
//    double lower_bound = 1e-10;
//    double upper_bound = 0.9;
//
//    lower_bound = 1e-10;
//    upper_bound = 0.9;
//
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<> rand_real(1,9);
//    std::uniform_int_distribution<> lower(6, 10);
//    std::uniform_int_distribution<> upper(1, 4);
//    lower_bound = rand_real(gen)*pow(10, -lower(gen));
//    upper_bound = rand_real(gen)*pow(10, -upper(gen));
//
//    lower_bound = 1e-10;
//    upper_bound = 0.9;
//
//    if (num_category == 2) {
//        parameters = {upper_bound, lower_bound};
//        cache_parameters = {upper_bound, lower_bound};
//    }
//    else {
//        std::cout << "Not yet implemented for more than 2 categories" << std::endl;
//        exit(222);
//        //TODO: Should throw exception instead of exit, this will do for now
//    }
//
////    em_model_ptr->at(0)->
//}


void EmAlgorithmMutation::InitialiseSummaryStat() {


    temp_stats = std::vector<std::vector<double>>(num_category);
    for (size_t i = 0; i < num_category; ++i) {
        all_em_stats.emplace_back(new EmSummaryStatMutation());
        temp_stats[i] = std::vector<double>(all_em_stats[i]->GetStatCount());
    }


}


void EmAlgorithmMutation::ExpectationStepCustom(size_t data_index, size_t category_index,
        double &sum_prob, std::vector<double> &temp_stat) {
    std::cout << "Error!! should NOT call ExpectationStepCustom here" << std::endl;
    exit(40);
}

