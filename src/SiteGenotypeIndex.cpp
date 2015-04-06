/*
 * site_genotypes_index.cc
 *
 *  Created on: 4/5/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/5/15.
//

#include "SiteGenotypeIndex.h"

/*
 * sequence_prob.cc
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#include <iostream>
#include <stdint.h>
#include <strings.h>
#include "site_genotypes.h"

const bool DEBUG = false;



SiteGenotypesIndex::~SiteGenotypesIndex() {
}

int SiteGenotypesIndex::GetDescendantCount() {
    return descendant_count;
}


void SiteGenotypesIndex::SetDescendantIndex(int des, int index) {
    descendant_index[des] = index;
}

int SiteGenotypesIndex::GetDescendantIndex(int des) {
    return descendant_index[des];

}

const std::vector<int>&SiteGenotypesIndex::GetDescendantIndex() {
    return descendant_index;
}

void SiteGenotypesIndex::SortIndex() {
//	descendant_index
//	std::algorithm::sort
    std::sort (descendant_index.begin(), descendant_index.end());

}
