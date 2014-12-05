#include <iostream>
#include "em_model_mutation.h"

EmModelMutation::EmModelMutation(EvolutionModel &evo_model0) : evo_model(&evo_model0){

    MutationRate rate = evo_model->GetMutationRate();
    cout << rate.prob << "\t" << rate.one_minus_p << endl;

}


void EmModelMutation::UpdateParameter(double param) {
    evo_model->UpdateMu(param);
//    MutationRate mutation_rate = evo_model->GetMutationRate();
//    cout << "IN EmModelMutation: updateing: " << param << "\t" << mutation_rate.prob << "\t" << mutation_rate.one_minus_p << endl;
}



EvolutionModel *EmModelMutation::GetEvoModel(){
    return evo_model;
}


void EmModelMutation::GetUpdatedInfo(){

    MutationMatrix transition_matrix_a_to_d = evo_model->GetTranstionMatirxAToD();
    MutationRate mutation_rate = evo_model->GetMutationRate();
    cout << "Called GetUpdateInfo with downcast: " << mutation_rate.prob << "\t" << mutation_rate.one_minus_p << endl;
}