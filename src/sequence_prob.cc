/*
 * sequence_prob.cc
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#include <iostream>
#include "sequence_prob.h"



const int ANCESTOR_COUNT = 10;
const int BASE_COUNT = 4;

Eigen::IOFormat nice_row(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");

const bool DEBUG = false;




SequenceProb::SequenceProb(ModelInput const site_data,  ModelParams const model_params) {

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

//	DiploidProbs num_genotypes = ancestor_genotypes;

    for (i = 1; i < site_data.all_reads.size(); ++i) {
        ReadData data = site_data.all_reads[i];
        all_descendant_data.push_back(data);

        HaploidProbs p = HaploidSequencing(data);
        all_descendant_genotypes.push_back(p);

    }

    descendant_count = site_data.all_reads.size() - 1;

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

DiploidProbs SequenceProb::DiploidSequencing(ReadData data) {
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

HaploidProbs SequenceProb::HaploidSequencing(ReadData data) {
	HaploidProbs result;
    
    double alphas_total = (1.0 - phi_haploid) / phi_haploid;
	for (int i : { 0, 1, 2, 3 }) {
		double alphas[4];
		for (int k : { 0, 1, 2, 3 }) {
			if (k == i)
				alphas[k] = (1.0 - error_prob) * alphas_total;
			else
				alphas[k] = error_prob / 3.0 * alphas_total;
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