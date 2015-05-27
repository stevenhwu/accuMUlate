#include <iostream>
#include "em_model_mutation_v1.h"


//EmModel(const EmModel &obj)
//{
//    cout << "Copy constructor allocating ptr." << endl;
//    ptr = new int;
//    *ptr = *obj.ptr; // copy the value
//}

EmModelMutationV1::EmModelMutationV1(const EmModelMutationV1 &em_model) {

    evo_model = em_model.evo_model->Clone().release();//get();
//    evo_model = em_model.evo_model->Clone2(); //CHECK: initial test no memory leak, new without delete?

    //MutationRate rate = evo_model->GetMutationRate();
    //std::cout << "Copy Constructor EmModelEvolutionV1: "<< rate.prob << "\t" << rate.one_minus_p << std::endl;
//TODO: Implemente rule of three,   Copy assignment operator!!

}


EmModelMutationV1::EmModelMutationV1(EvolutionModel &evo_model0) : evo_model(&evo_model0){

    //MutationRate rate = evo_model->GetMutationRate();
    //std::cout << "Constructor EmModelEvolutionV1: " << rate.prob << "\t" << rate.one_minus_p << std::endl;

}


void EmModelMutationV1::UpdateParameter(double param) {

//    cout << "IN EmModelMutationV1: updateing: " << param << "\t" << mutation_rate.prob << "\t" << mutation_rate.one_minus_p << endl;
    evo_model->UpdateOneMinusExpBeta(param);


}


void EmModelMutationV1::GetParameterInfo(){

    MutationMatrix transition_matrix_a_to_d = evo_model->GetTranstionMatirxAToD();
    double mutation_rate = evo_model->GetMutationRate();
    MutationProb mutation_prob = evo_model->GetMutationProb();
    std::cout << "Called GetUpdateInfo: " << mutation_rate << "\t" << (1.0-mutation_rate) <<  "\t" << mutation_prob.GetExpBeta()<< std::endl;
}


EvolutionModel *EmModelMutationV1::GetEvoModel() {
    return evo_model;
}

void EmModelMutationV1::UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat, double &log_likelihood_scaler) {
    std::cout << "Error!! Not yet implemented" << std::endl;
    std::exit(44);
}

size_t EmModelMutationV1::GetDataCount() {
    std::cout << "Error!! Not yet implemented" << std::endl;
    std::exit(44);
}

//EvolutionModel * EmModelMutationV1::GetEvoModel() const {
//    return evo_model;
//}
//
//std::unique_ptr<EvolutionModel> EmModelMutationV1::CopyEvoModel()const {
//    std::unique_ptr<EvolutionModel> e = evo_model->Clone();
//    return e;
//}
//
//EvolutionModel* EmModelMutationV1::CopyEvoModel2()const {
//    EvolutionModel* e = evo_model->Clone2();
//    return e;
//}
//
//EmModelMutationV1* EmModelMutationV1::GetModel() {
//    return this;
//}
