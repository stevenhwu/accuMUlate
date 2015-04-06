#include "em_algorithm_mutation.h"
#include "em_summary_stat_mutation.h"


EmAlgorithmMutation::EmAlgorithmMutation(std::vector<std::unique_ptr<EmModel>> &model_ptr) : EmAlgorithm(model_ptr) {

    InitWithModel();

}


EmAlgorithmMutation::~EmAlgorithmMutation() {

}

void EmAlgorithmMutation::Run() {

//    em_stat_local_single->Print();
    size_t i = 0;
    bool isConverged = true;
    while(isConverged){
        ExpectationStepModelPtr();
        MaximizationStep();
        isConverged = EmStoppingCriteria(i);
        i++;
    }
}


void EmAlgorithmMutation::InitialiseParameters() {
    double lower_bound = 1e-10;
    double upper_bound = 0.9;

    if (num_category == 2) {
        parameters = {upper_bound, lower_bound};
    }
    else {
        std::cout << "Not yet implemented for more than 2 categories" << std::endl;
        exit(222);
        //TODO: Should throw exception instead of exit, this will do for now
    }

//    em_model_ptr->at(0)->
}


void EmAlgorithmMutation::InitialiseSummaryStat() {

//    em_stat_local_single = std::unique_ptr<EmSummaryStat>(new EmSummaryStatMutation());
//    em_stat_local_single->Print();

    temp_stats = std::vector<std::vector<double>>(num_category);
    for (size_t i = 0; i < num_category; ++i) {
        all_em_stats.emplace_back(new EmSummaryStatMutation());
//        all_em_stats.emplace_back(new EmSummaryStatMutation());
        temp_stats[i] = std::vector<double>(all_em_stats[i]->GetStatCount());
    }

//    std::cout << "temp_stat_size: " << temp_stats.size() << std::endl;
//    for (int i = 0; i < temp_stats.size(); ++i) {
//        std::cout << "temp_stat_size[i]: " << temp_stats[i].size() << std::endl;
//        for (int j = 0; j < temp_stats[i].size(); ++j) {
//            std::cout << temp_stats[i][j] << "\t";
//        }
//        std::cout << std::endl;
//    }


}


void EmAlgorithmMutation::ExpectationStepCustom(size_t data_index, size_t category_index,
        double &sum_prob, std::vector<double> &temp_stat) {
    std::cout << "Error!! should NOT call ExpectationStepCustom here" << std::endl;
    exit(40);
}
