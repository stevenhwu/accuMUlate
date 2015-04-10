/*
 * sequence_prob.cc
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#include <iostream>
#include <stdint.h>
#include <strings.h>
#include "sequence_prob_v1.h"



const int ANCESTOR_COUNT = 10;
const int BASE_COUNT = 4;

//Eigen::IOFormat nice_row(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");

const bool DEBUG = false;



SequenceProb::SequenceProb(ModelInput const &site_data,  ModelParams const &model_params) {

    phi_haploid = model_params.phi_haploid;
    phi_diploid = model_params.phi_diploid;
    error_prob = model_params.error_prob;
    theta = model_params.theta;
	frequency_prior = model_params.nuc_freq;

    size_t i=0;
    DiploidProbs pop_genotypes = DiploidPopulation(site_data.reference);
    ancestor_data = site_data.all_reads[i];
	ancestor_genotypes = DiploidSequencing(ancestor_data);
	ancestor_genotypes *= pop_genotypes;

	descendant_count = site_data.all_reads.size() - 1;
	descendant_index.assign(descendant_count, -1);

	all_descendant_data.reserve(descendant_count);
	all_descendant_genotypes.reserve(descendant_count);


    for (i = 1; i < (descendant_count+1); ++i) {
        ReadData data = site_data.all_reads[i];
        all_descendant_data.push_back(data);
//		all_descendant_data[i-1] =data ;

//        HaploidProbs p = HaploidSequencing(data);
		all_descendant_genotypes.emplace_back(HaploidSequencing(data));
//        all_descendant_genotypes.push_back(p);
//		all_descendant_genotypes[i-1]=p;

    }

//	for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
//		int index16 = LookupTable::index_converter_10_to_16[index10];
//
//		double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[index10] *  prod_prob_ancestor;
//		all_sequence_prob[site_index].GetAncestorGenotypesMultiplyPrior(index10);
//		prob_reads += prob_reads_given_a;


////	std::cout << "==================="<< std::endl;
////
//	for (auto item : temp_map) {
//
//		condense_genotype.emplace_back(std::make_pair(item.first, item.second));
////		auto a = std::make_pair(item.first, item.second);
////		std::cout << item.first << "\t" << item.second << "\t" << a.first << "\t" << a.second <<std::endl;
//	}
////
////	std::cout << "New count:\t" << condense_genotype.size() << std::endl;
////	for (auto item : condense_genotype) {
////		std::cout << item.first << "\t" << item.second <<std::endl;
////	}


}


SequenceProb::~SequenceProb() {
}




DiploidProbs SequenceProb::DiploidPopulation(int ref_allele) { //TODO: refactor out, noly need to do this 4 times
	ReadData d;
	DiploidProbs result;

	double alphas[4];
	for(int i : {0,1,2,3}){

        alphas[i] = theta * frequency_prior[i];
//		printf("A:%f\n",alphas[i]);
	}
//	printf("ref:%d\n", ref_allele);
	for(int i : {0,1,2,3}) {
		for(int j=0;j<i;++j) {
			d.key = 0;
			d.reads[ref_allele] = 1;
			d.reads[i] += 1;
			d.reads[j] += 1;
//			printf("  %d %d: %u %u %u %u\n", i, j, d.reads[0], d.reads[1], d.reads[2], d.reads[3]);
			result[i+j*4] = DirichletMultinomialLogProbability(alphas, d);
			result[j+i*4] = result[i+j*4];
		}
		d.key = 0;
		d.reads[ref_allele] = 1;
		d.reads[i] += 2;
//		printf("D:%d %d: %u %u %u %u\n", i, i, d.reads[0], d.reads[1], d.reads[2], d.reads[3]);
		result[i+i*4] = DirichletMultinomialLogProbability(alphas, d);
	}
//	std::cout << result << std::endl;
	return result.exp();
}

DiploidProbs SequenceProb::DiploidSequencing(ReadData const &data) {
	DiploidProbs result;
    double alphas_total = (1.0 - phi_diploid) / phi_diploid;
	for (int i : { 0, 1, 2, 3 }) {
        for (int j = 0; j < i; ++j) {
			double alphas[4];
			for (int k : { 0, 1, 2, 3 }) {
				if (k == i || k == j)
					alphas[k] = (0.5 - error_prob / 3.0) * alphas_total;
				else
					alphas[k] = (error_prob / 3.0) * alphas_total;
			}
			result[i * 4 + j] = DirichletMultinomialLogProbability(alphas, data);
			result[j * 4 + i] = result[i * 4 + j];
		}
		double alphas[4];
		for (int k : { 0, 1, 2, 3 }) {
			if (k == i)
				alphas[k] = (1.0 - error_prob) * alphas_total;
			else
				alphas[k] = error_prob / 3.0 * alphas_total;
		}
		result[i * 4 + i] = DirichletMultinomialLogProbability(alphas, data);
	}


	double scale = result.maxCoeff();
	return (result - scale).exp();
//    result -= scale;
//    result = ScaleLogArray(result);
//    return result.exp();
}

HaploidProbs SequenceProb::HaploidSequencing(ReadData const &data) {
	HaploidProbs result;
    
    double alphas_total = (1.0 - phi_haploid) / phi_haploid;
	for (int i : { 0, 1, 2, 3 }) {
		double alphas[4];
		std::fill(alphas, alphas+4, error_prob / 3.0 * alphas_total);
		for (int k : { 0, 1, 2, 3 }) {
			if (k == i) {
				alphas[k] = (1.0 - error_prob) * alphas_total;
				break;
			}
//			else[k] = error_prob / 3.0 * alphas_total;
		}
		result[i] = DirichletMultinomialLogProbability(alphas, data);
	}

	double scale = result.maxCoeff();
	return (result - scale).exp();
//    result -= scale;
//    result = ScaleLogArray(result);
//    return result.exp();
}

// FIXME : make sure it only take Eigen::Array class, pass by value or ref here?? Turn out it's not the C++ way, need to use static_cast
template <class T>
T SequenceProb::ScaleLogArray(T result) {
    double scale = result.maxCoeff();
//    double normalise = result.sum();
//    result /=  normalise;

	result = (result - scale).exp();
    result = result.exp();

    return result;
}

std::vector<HaploidProbs> SequenceProb::GetDescendantGenotypes() {
    return all_descendant_genotypes;
}

HaploidProbs SequenceProb::GetDescendantGenotypes(int descent_index) {

    if(descent_index >= descendant_count || descent_index < 0){
        //TODO: throw error
    }

    return all_descendant_genotypes[descent_index];
}

DiploidProbs SequenceProb::GetAncestorGenotypes() {
    return ancestor_genotypes;
}


//HaploidProbs SequenceProb::HaploidProbsFactory(ReadData const &data) {
//	return SequenceProb::HaploidSequencing(data);
//}

std::array<DiploidProbs, 4> SequenceProb::DiploidPopulationFactory(ModelParams const model_params) {


	double theta = model_params.theta;
	std::vector<double> frequency_prior = model_params.nuc_freq;

	double alphas[4];
	for(int i : {0,1,2,3}){
		alphas[i] = theta * frequency_prior.at(i);
	}

	ReadData d;
	DiploidProbs result;
	std::array<DiploidProbs, 4> cache_results;
	for (int ref = 0; ref < 4; ++ref) {
		for (int i : {0, 1, 2, 3}) {
			for (int j = 0; j < i; ++j) {
				d.key = 0;
				d.reads[ref] = 1;
				d.reads[i] += 1;
				d.reads[j] += 1;
				result[i + j * 4] = DirichletMultinomialLogProbability(alphas, d);
				result[j + i * 4] = result[i + j * 4];
			}
			d.key = 0;
			d.reads[ref] = 1;
			d.reads[i] += 2;
			result[i + i * 4] = DirichletMultinomialLogProbability(alphas, d);
		}
		cache_results[ref] = result.exp();
	}

	return cache_results;
}

int SequenceProb::GetDescendantCount() {
	return descendant_count;
}

void SequenceProb::printReadData(ReadData read_data) {
	printf("%d %d %d %d\n", read_data.reads[0], read_data.reads[1], read_data.reads[2], read_data.reads[3]);

}

void SequenceProb::PrintReads(ReadData data) {
	printf("%d %d %d %d\n", data.reads[0], data.reads[1], data.reads[2], data.reads[3]);

}

ReadData SequenceProb::GetDescendantReadData(int descent_index) {
	return all_descendant_data[descent_index];
}
ReadDataVector SequenceProb::GetDescendantReadData() {
	return all_descendant_data;
}

ReadDataVector const & SequenceProb::GetDescendantReadData2() {
	return all_descendant_data;
}

ReadDataVector const * SequenceProb::GetDescendantReadData3() {
	return &all_descendant_data;
}

void SequenceProb::GetDescendantReadDataCOPY(ReadDataVector &all_descendant_data2) {
	all_descendant_data2 = all_descendant_data;
}

uint64_t SequenceProb::GetDescendantReadDataKey(int descent_index) {
	return all_descendant_data[descent_index].key;
}

void SequenceProb::SetDescendantIndex(int des, int index) {
	descendant_index[des] = index;
}

int SequenceProb::GetDescendantIndex(int des) {
	return descendant_index[des];

}

const std::vector<int>& SequenceProb::GetDescendantIndex() {
	return descendant_index;
}

void SequenceProb::SortIndex() {
//	descendant_index
//	std::algorithm::sort
	std::sort (descendant_index.begin(), descendant_index.end());

}

SequenceProb::SequenceProb() {

}
