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



//std::vector<std::unique_ptr<EmData>>
EmAlgorithmMutationV1::EmAlgorithmMutationV1(int num_category0, std::vector<SiteProb> &em_data0, EvolutionModel &em_model0,
        std::vector<std::unique_ptr<EmData>> &d_ptr, EmModel &m)
        : EmAlgorithm(num_category0, d_ptr, m), em_data_old(em_data0), em_model_old(&em_model0)
//          EmAlgorithm::em_data_ptr(&d_ptr), em_model(&m) {
{

    if (num_category != 2) {
        std::cout << "Not yet implemented for more than 2 categories: input_cat: " << "\t" << num_category<< std::endl;
        exit(222);
    }

    InitWithData();


    MutationProb mutation_prob = em_model_old->GetMutationProb();//FIXME: What to do here?? maybe in EmSummayrStat??

    double lower_bound = mutation_prob.ConvertExpBetaToMu(1e-12);
    double upper_bound = mutation_prob.ConvertExpBetaToMu(1e-1);
    all_probs_old = Eigen::ArrayXXd::Zero(num_category, site_count);
    parameters_old = std::vector<double>(num_category);
    if (num_category == 2) {
        parameters_old = {upper_bound, lower_bound};
    }
    all_stats_same = std::vector<double>(num_category, 0);
    all_stats_diff = std::vector<double>(num_category, 0);


    std::cout << "=========== Done Constructor smart pointer\n";

}


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
//    max_ite_count = 3;
//
////    auto site2 = em_data[0];
////    EmData *site3 = em_data[0];
////
////    site2->Test(4);
////    site3->Test(5);
//
//
//    InitWithData();
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
//    max_ite_count = 10;
//
//    InitWithData();
//}


EmAlgorithmMutationV1::~EmAlgorithmMutationV1() {

}

void EmAlgorithmMutationV1::Run() {
em_stat_local_single->print();
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

    em_stat_local_single->print();
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


void EmAlgorithmMutationV1::InitialiseSummaryStat() {

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




