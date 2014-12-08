

#include <iostream>
#include "em_data_mutation.h"
#include "em_model_mutation.h"


EmDataMutation::EmDataMutation(SequenceProb &sequence_prob, EvolutionModel &evo_model){
    site = SiteProb(sequence_prob, evo_model);
}


EmDataMutation::~EmDataMutation() {}


void EmDataMutation::UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat) {


    summaryStat->print();
    prob = 0;
    double stat_same = 0;
    double stat_diff = 0;
    //relpace with
//    EmSummaryStatMutation local;


    site.CalculateAncestorToDescendant(prob, stat_same, stat_diff);

    summaryStat->SetStats(vector<double> {stat_same, stat_diff});

//    std::cout << "call EMDATA2: " << stat_same << "\t" << stat_diff << "\t" << prob <<std::endl;
//    site->UpdateSummaryStat(sum_prob, local);
//    site.CalculateAncestorToDescendant(sum_prob, stat_same, stat_diff);
}


void EmDataMutation::UpdateEmModel(EmModel *em_model) {

//    cout << "In EmDataMutation, UpdateModel: " << endl;
//    EmModelMutation *em = dynamic_cast<EmModelMutation*> (em_model);
    EmModelMutation *em_model_mutation = static_cast<EmModelMutation*> (em_model);
    EvolutionModel *evo_model = em_model_mutation->GetEvoModel();
    site.UpdateModel(*evo_model);
    //TODO: redo static_cast with visitor??? ?typeid?

}

void EmDataMutation::UpdateEmModel(std::unique_ptr<EmModelMutation> &em_model) {

//    default_delete<EmModel> &deleter = em_model.get_deleter();
//    std::unique_ptr<EmDataMutation> em_model_mutation = static_unique_ptr_cast<EmModelMutation, EmModel, default_delete<EmModel> >(em_model);
    EvolutionModel *evo_model = em_model->GetEvoModel();
    site.UpdateModel(*evo_model);
}



/*
FIXME: This is getting stupid, fix it later
template<typename Derived, typename Base, typename Del>
std::unique_ptr<Derived, Del>
static_unique_ptr_cast( std::unique_ptr<Base, Del>&& p )
{
    auto d = static_cast<Derived *>(p.release());
    return std::unique_ptr<Derived, Del>(d, std::move(p.get_deleter()));
}

template<typename Derived, typename Base, typename Del>
std::unique_ptr<Derived, Del>
dynamic_unique_ptr_cast( std::unique_ptr<Base, Del>&& p )
{
    if(Derived *result = dynamic_cast<Derived *>(p.get())) {
        p.release();
        return std::unique_ptr<Derived, Del>(result, std::move(p.get_deleter()));
    }
    return std::unique_ptr<Derived, Del>(nullptr, p.get_deleter());
}


void EmDataMutation::UpdateEmModel(std::unique_ptr<EmModel> &em_model) {

    default_delete<EmModel> &deleter = em_model.get_deleter();
    std::unique_ptr<EmDataMutation> em_model_mutation = static_unique_ptr_cast<EmModelMutation, EmModel, default_delete<EmModel> >(em_model);

}
*/
void EmDataMutation::UpdateEmModel(std::unique_ptr<EmModelBinomial> &em_model) {

}
