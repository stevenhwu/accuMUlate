/*
 * site_genotypes_index.cc
 *
 *  Created on: 4/5/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/5/15.
//

#ifndef _SITE_GENOTYPES_INDEX_H_
#define _SITE_GENOTYPES_INDEX_H_

#include <vector>
#include <iostream>
#include "model.h"
#include "mutation_prob.h"
#include "site_prob.h"
#include "site_genotypes.h"
#include "distributions/DirichletMultinomialDistribution.h"
#include "evolution_models/EvolutionModel.h"
#include "lookup.h"
#include "sequence_prob_v1.h"
#include "constant.h"

class SiteGenotypesIndex {
public:


    SiteGenotypesIndex(){}

    SiteGenotypesIndex(uint16_t ref, int descendant_count0): reference(ref), descendant_count(descendant_count0){
//        all_descendant_genotypes.reserve(descendant_count);
        descendant_index.assign(descendant_count, -1);
    }


    uint16_t GetReference(){
        return reference;
    }

    void SetAncestorIndex(int index){
        ancestor_index = index;
    }
    int GetAncestorIndex(){
        return ancestor_index;
    }

    ~SiteGenotypesIndex();


    int GetDescendantCount();


    void SetDescendantIndex(int des, int index);
    int GetDescendantIndex(int des);


    const std::vector<int>& GetDescendantIndex();

    void SortIndex();




private:


    size_t descendant_count;
    uint16_t reference;
    std::vector<int> descendant_index;
    int ancestor_index;

};


#endif //_SITE_GENOTYPES_INDEX_H_
