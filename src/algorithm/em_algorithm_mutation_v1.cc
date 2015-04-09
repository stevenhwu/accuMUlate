#include "em_algorithm_mutation_v1.h"

//data_count: 41321 ite:33  Total EM Time: 296.134882

EmAlgorithmMutationV1::EmAlgorithmMutationV1(int num_category0,
        std::vector<std::unique_ptr<EmData>> &d_ptr, EmModelMutationV1 &m)
        : EmAlgorithm(num_category0, d_ptr, m)
//          EmAlgorithm::em_data_ptr(&d_ptr), em_model(&m) {
{

    for (size_t i = 0; i < num_category; ++i) {
        em_model.emplace_back(new EmModelMutationV1(m));
    }


    std::cout << "data_count: " << em_data_ptr->size() << std::endl;
    std::cout << "=========== Done Constructor with smart pointer and copy constructor\n";
    InitWithData();


}



EmAlgorithmMutationV1::EmAlgorithmMutationV1(
        std::vector<std::unique_ptr<EmData>> &d_ptr, std::vector<std::unique_ptr<EmModel>> &m)
        : EmAlgorithm(d_ptr, m)
//          EmAlgorithm::em_data_ptr(&d_ptr), em_model(&m) {
{


    InitWithData();
    std::cout << "=========== Done Constructor smart pointer x 2\n";

}




EmAlgorithmMutationV1::~EmAlgorithmMutationV1() {

}

void EmAlgorithmMutationV1::RunEM() {
//    em_stat_local_single->Print();
    size_t i = 0;
    bool isConverged = true;
    while(isConverged){
        ExpectationStepModel();
        MaximizationStep();
        isConverged = EmStoppingCriteria(i);
        i++;
    }

}

void EmAlgorithmMutationV1::ExpectationStepCustom(size_t data_index, size_t category_index,
        double &sum_prob, std::vector<double> &temp_stat) {

    em_data_ptr->at(data_index)->UpdateEmModel( em_model[category_index].get() );
    em_data_ptr->at(data_index)->UpdateSummaryStat(sum_prob, temp_stat);

}



void EmAlgorithmMutationV1::Run2() {


    size_t i = 0;
    bool isConverged = true;
    while (isConverged) {
        ExpectationStepModelPtr();
        MaximizationStep();
        isConverged = EmStoppingCriteria(i);
        i++;
    }
}

void EmAlgorithmMutationV1::InitialiseParameters() {
    double lower_bound = 1e-1;
    double upper_bound = 0.5;

    lower_bound = 1e-10;
    upper_bound = 0.9;

    if (num_category == 2) {
        parameters = {upper_bound, lower_bound};
    }
    else {
        std::cout << "Not yet implemented for more than 2 categories" << std::endl;
        exit(222);
        //TODO: Should throw exception instead of exit, this will do for now
    }


}


void EmAlgorithmMutationV1::InitialiseSummaryStat() {

//    em_stat_local_single = std::unique_ptr<EmSummaryStat>(new EmSummaryStatMutation());
//    em_stat_local_single->Print();

    temp_stats = std::vector<std::vector<double>>(num_category);
    for (size_t i = 0; i < num_category; ++i) {
        all_em_stats.emplace_back(new EmSummaryStatMutation());
//        all_em_stats.emplace_back(new EmSummaryStatMutation());
        temp_stats[i] = std::vector<double>(all_em_stats[i]->GetStatCount());
    }


//    for (size_t i = 0; i < num_category; ++i) {
//        EmSummaryStatMutation a;
//        all_em_stats_nonptr.push_back(a);
////        all_em_stats_nonptr[i]->SetStats(std::vector<double> {i+1.0, i+1.0});
////        all_em_stats_nonptr[i]->Print();
//    }
//    for (size_t i = 0; i < num_category; ++i) {
//
////        all_em_stats_nonptr[i]->Print();
//    }

}




