

#include <iostream>
#include "em_data_mutation.h"


EmDataMutation::EmDataMutation(SequenceProb &sequence_prob, EvolutionModel &evo_model){
    site = SiteProb(sequence_prob, evo_model);
}


EmDataMutation::~EmDataMutation() {}





void EmDataMutation::UpdateLikelihood(double prob, EmSummaryStat &summaryStat) {


    std::cout << "call EMDATA2" << std::endl;
    summaryStat.print();
}




