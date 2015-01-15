#include "em_algorithm_mutation.h"
#include "em_summary_stat_mutation_v1.h"


//std::vector<std::unique_ptr<EmData>>
//EmAlgorithmMutation::EmAlgorithmMutation(int num_category0,
//        std::vector<std::unique_ptr<EmData>> &d_ptr, EmModelMutationV1 &m)
//        : EmAlgorithm(num_category0, d_ptr, m)
////          EmAlgorithm::em_data_ptr(&d_ptr), em_model(&m) {
//{
//
//    if (num_category != 2) {
//        std::cout << "Not yet implemented for more than 2 categories: input_cat: " << "\t" << num_category<< std::endl;
//        exit(222);
//    }
//
//    for (size_t i = 0; i < num_category; ++i) {
//        em_model.emplace_back(new EmModelMutationV1(m));
//    }
//
//    em_count = 1000;
//    std::cout << "data_count: " << em_data_ptr->size() << std::endl;
//    Init();
//
////    em_model_ptr->at(0)->GetParameterInfo();
////    em_model_ptr->at(1)->GetParameterInfo();
////    em_model[0]->UpdateParameter(parameters[0]);
////    em_model[1]->UpdateParameter(parameters[1]);
//
//
//}

EmAlgorithmMutation::EmAlgorithmMutation(std::vector<std::unique_ptr<EmModel>> &model_ptr) : EmAlgorithm(model_ptr) {


    InitModelOnly();

}


//
//EmAlgorithmMutation::EmAlgorithmMutation(
//        std::vector<std::unique_ptr<EmData>> &d_ptr, std::vector<std::unique_ptr<EmModel>> &m)
//        : EmAlgorithm(d_ptr, m)
////          EmAlgorithm::em_data_ptr(&d_ptr), em_model(&m) {
//{
//
//    if (num_category != 2) {
//        std::cout << "Not yet implemented for more than 2 categories: input_cat: " << "\t" << num_category<< std::endl;
//        exit(222);
//    }
//    em_count = 10;
//
//    Init();
//
//
//    std::cout << "=========== Done Constructor smart pointer x 2\n";
//
//}
//
//
//
////std::vector<std::unique_ptr<EmData>>
//EmAlgorithmMutation::EmAlgorithmMutation(int num_category0, std::vector<SiteProb> &em_data0, EvolutionModel &em_model0,
//        std::vector<std::unique_ptr<EmData>> &d_ptr, EmModel &m)
//        : EmAlgorithm(num_category0, d_ptr, m), em_data_old(em_data0), em_model_old(&em_model0)
////          EmAlgorithm::em_data_ptr(&d_ptr), em_model(&m) {
//{
////    cout << "\n============\nEM  half way Constructor\n";
//
////    em_model->UpdateParameter(2);
//    if (num_category != 2) {
//        std::cout << "Not yet implemented for more than 2 categories: input_cat: " << "\t" << num_category<< std::endl;
//        exit(222);
//    }
//    em_count = 10;
//
//    Init();
//
//
//    MutationProb mutation_prob = em_model_old->GetMutationProb();//FIXME: What to do here?? maybe in EmSummayrStat??
//
//    double lower_bound = mutation_prob.ConvertExpBetaToMu(1e-12);
//    double upper_bound = mutation_prob.ConvertExpBetaToMu(1e-1);
//    all_probs_old = Eigen::ArrayXXd::Zero(num_category, site_count);
//    parameters_old = std::vector<double>(num_category);
//    if (num_category == 2) {
//        parameters_old = {upper_bound, lower_bound};
//    }
//    all_stats_same = std::vector<double>(num_category, 0);
//    all_stats_diff = std::vector<double>(num_category, 0);
//
//
//    std::cout << "=========== Done Constructor smart pointer\n";
//
//}


//
//EM::EM(int num_category0, vector<SiteProb> &em_data0, EvolutionModel &em_model0, vector<EmData *> &d, EmModel &m)
//        : num_category(num_category0), em_data_old(em_data0), em_model_old(&em_model0), em_data(d), em_model(&m) {
//
////    cout << "\n============\nEM  half way Constructor\n";
//
////    em_model->UpdateParameter(2);
//    if (num_category != 2) {
//        cout << "Not yet implemented for more than 2 categories" << endl;
//        exit(222);
//    }
//
//    site_count = em_data.size();
//    num_category = num_category;
//    em_count = 3;
//
////    auto site2 = em_data[0];
////    EmData *site3 = em_data[0];
////
////    site2->Test(4);
////    site3->Test(5);
//
//
//    Init();
//    cout << "===========\nDone Constructor raw pointer\n";
//}


//EM::EM(int num_category0, vector<EmData> em_data0, EmModel &em_model0) : num_category(num_category0), em_data(em_data0), em_model(&em_model0){
//
//    cout << "EM Constructor\n";
//    if (num_category != 2) {
//        cout << "Not yet implemented for more than 2 categories" << endl;
//        exit(222);
//    }
//
//    site_count = em_data_old.size();
//    num_category = num_category;
//    em_count = 10;
//
//    Init();
//}


EmAlgorithmMutation::~EmAlgorithmMutation() {

}

void EmAlgorithmMutation::Run() {
em_stat_local_single->print();
    size_t i = 0;
    bool isConverged = true;
    while(isConverged & (i< em_count) ){
//        cout << "EM ite: " << i << endl;
        ExpectationStep();

//        oldEStep();
//        for (int r = 0; r < num_category; ++r) {
//            if (all_em_stats[r]->GetStat(0) != all_stats_same[r] || all_em_stats[r]->GetStat(1) != all_stats_diff[r]) {
//                cout << "DIFF!!! inner " << all_em_stats[r]->GetStat(0) << "\t" << all_stats_same[r] << "\t" << all_em_stats[r]->GetStat(1) << "\t" << all_stats_diff[r] << endl;
//                exit(99);
//            }
//        }


        MaximizationStep();


//        oldMStep();
        isConverged = EmStoppingCriteria(i);
        i++;
    }

}



void EmAlgorithmMutation::Run2() {

    for (size_t i = 0; i < 100; ++i) {
        std::cout << "EM2 ite: " << i << std::endl;
        ExpectationStep2();
        MaximizationStep();
    }
}
void EmAlgorithmMutation::ExpectationStepCustom(size_t data_index, size_t category_index,
        double &sum_prob, std::vector<double> &temp_stat) {

    em_data_ptr->at(data_index)->UpdateEmModel( em_model[category_index].get() );
    em_data_ptr->at(data_index)->UpdateSummaryStat(sum_prob, temp_stat);


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


}


void EmAlgorithmMutation::InitialiseSummaryStat() {

    em_stat_local_single = std::unique_ptr<EmSummaryStat>(new EmSummaryStatMutationV1());
    em_stat_local_single->print();

    temp_stats = std::vector<std::vector<double>>(num_category);
    for (size_t i = 0; i < num_category; ++i) {
        all_em_stats.emplace_back(new EmSummaryStatMutationV1());
//        all_em_stats.emplace_back(new EmSummaryStatMutationV1());
        temp_stats[i] = std::vector<double>(em_stat_local_single->GetStatCount());
    }


//    for (size_t i = 0; i < num_category; ++i) {
//        EmSummaryStatMutationV1 a;
//        all_em_stats_nonptr.push_back(a);
////        all_em_stats_nonptr[i]->SetStats(std::vector<double> {i+1.0, i+1.0});
////        all_em_stats_nonptr[i]->print();
//    }
//    for (size_t i = 0; i < num_category; ++i) {
//
////        all_em_stats_nonptr[i]->print();
//    }

}


