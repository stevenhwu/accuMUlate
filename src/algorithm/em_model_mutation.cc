#include <iostream>
#include "em_model_mutation.h"


//EmModel(const EmModel &obj)
//{
//    cout << "Copy constructor allocating ptr." << endl;
//    ptr = new int;
//    *ptr = *obj.ptr; // copy the value
//}

EmModelMutation::EmModelMutation(const EmModelMutation &em_model) {

    evo_model = em_model.evo_model->Clone().release();//get();
//    evo_model = em_model.evo_model->Clone2(); //CHECK: initial test no memory leak, new without delete?

    MutationRate rate = evo_model->GetMutationRate();
    cout << "Copy Constructor EmModelEvolution: "<< rate.prob << "\t" << rate.one_minus_p << endl;
//FIXME: Implemente rule of three,   Copy assignment operator!!

}


EmModelMutation::EmModelMutation(EvolutionModel &evo_model0) : evo_model(&evo_model0){

    MutationRate rate = evo_model->GetMutationRate();
    cout << rate.prob << "\t" << rate.one_minus_p << endl;

}


void EmModelMutation::UpdateParameter(double param) {

//    MutationRate mutation_rate = evo_model->GetMutationRate();
//    cout << "IN EmModelMutation: updateing: " << param << "\t" << mutation_rate.prob << "\t" << mutation_rate.one_minus_p << endl;
    evo_model->UpdateExpBeta(param);


}


void EmModelMutation::GetParameterInfo(){

    MutationMatrix transition_matrix_a_to_d = evo_model->GetTranstionMatirxAToD();
    MutationRate mutation_rate = evo_model->GetMutationRate();
    MutationProb mutation_prob = evo_model->GetMutationProb();
    cout << "Called GetUpdateInfo: " << mutation_rate.prob << "\t" << mutation_rate.one_minus_p <<  "\t" << mutation_prob.GetExpBeta()<< endl;
}


EvolutionModel * EmModelMutation::GetEvoModel() {
    return evo_model;
}
//EvolutionModel * EmModelMutation::GetEvoModel() const {
//    return evo_model;
//}
//
//std::unique_ptr<EvolutionModel> EmModelMutation::CopyEvoModel()const {
//    std::unique_ptr<EvolutionModel> e = evo_model->Clone();
//    return e;
//}
//
//EvolutionModel* EmModelMutation::CopyEvoModel2()const {
//    EvolutionModel* e = evo_model->Clone2();
//    return e;
//}
//
//EmModelMutation* EmModelMutation::GetModel() {
//    return this;
//}
