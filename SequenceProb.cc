/*
 * SequencingProb.cc
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#include "SequenceProb.h"

#include <iostream>




//
//double TetMAProbability(const ModelParams &params, const ModelInput site_data) {
//	MutationMatrix m = MutationAccumulation(params, false);
//	MutationMatrix mt = MutationAccumulation(params, true);
//
//	MutationMatrix mn = m - mt;
//	DiploidProbs pop_genotypes = DiploidPopulation(params, site_data.reference);
//
//	auto it = site_data.all_reads.begin();
//	DiploidProbs anc_genotypes = DiploidSequencing(params, site_data.reference,
//			*it);
//	anc_genotypes *= pop_genotypes;
//	DiploidProbs num_genotypes = anc_genotypes;
//	for (++it; it != site_data.all_reads.end(); ++it) {
//		HaploidProbs p = HaploidSequencing(params, site_data.reference, *it);
//
//		anc_genotypes *= (m.matrix() * p.matrix()).array();
//		num_genotypes *= (mn.matrix() * p.matrix()).array();
//		exit(-1);
//	}
//
//	//cerr << "\n" << anc_genotypes/anc_genotypes.sum() << endl;
//
//	return 1.0 - num_genotypes.sum() / anc_genotypes.sum();
//}

SequenceProb::SequenceProb(const ModelParams &model_params,
		const ModelInput site_data, MutationProb muProb) {

	params = model_params; //TODO check copy method/ref

	UpdateMuProb(muProb);
	data = site_data;
	pop_genotypes = DiploidPopulation(site_data.reference);

	auto it = site_data.all_reads.begin();
	anc_genotypes = DiploidSequencing(*it);
	anc_genotypes *= pop_genotypes;
	DiploidProbs num_genotypes = anc_genotypes;
	for (++it; it != site_data.all_reads.end(); ++it) {
		HaploidProbs p = HaploidSequencing(*it);
		allHap.push_back(p);
	}
	UpdateLikelihood();

}

SequenceProb::~SequenceProb() {
}

void SequenceProb::UpdateMuProb(MutationProb muProb){

	mutation = muProb.GetMutation();
	non_mutation = muProb.GetNonMutation();
	UpdateLikelihood();
}

DiploidProbs SequenceProb::DiploidPopulation(int ref_allele) {
	ReadData d;
	DiploidProbs result;

	double alphas[4];
	for(int i : {0,1,2,3}){
		alphas[i] = params.theta*params.nuc_freq[i];
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
//			result[i+j*4] = 5;
//			result[j+i*4] = 3;
		}
		d.key = 0;
		d.reads[ref_allele] = 1;
		d.reads[i] += 2;
//		printf("D:%d %d: %u %u %u %u\n", i, i, d.reads[0], d.reads[1], d.reads[2], d.reads[3]);
		result[i+i*4] = DirichletMultinomialLogProbability(alphas, d);
	}
//	std::cout << result << std::endl;
//	Eigen::Matrix4d m;
//	for (int i = 0; i < result.rows(); ++i) {
//		m(i)= result[i];
//	}
//	std::cout << m<< std::endl;
	return result.exp();
}
DiploidProbs SequenceProb::DiploidSequencing(ReadData data) {
	DiploidProbs result;
	double alphas_total = (1.0 - params.phi_diploid) / params.phi_diploid;
	for (int i : { 0, 1, 2, 3 }) {
		for (int j = 0; j < i; ++j) {
			double alphas[4];
			for (int k : { 0, 1, 2, 3 }) {
				if (k == i || k == j)
					alphas[k] = (0.5 - params.error_prob / 3.0) * alphas_total;
				else
					alphas[k] = (params.error_prob / 3.0) * alphas_total;
			}
			result[i * 4 + j] = DirichletMultinomialLogProbability(alphas,
					data);
			result[j * 4 + i] = result[i * 4 + j];
		}
		double alphas[4];
		for (int k : { 0, 1, 2, 3 }) {
			if (k == i)
				alphas[k] = (1.0 - params.error_prob) * alphas_total;
			else
				alphas[k] = params.error_prob / 3.0 * alphas_total;
		}
		result[i * 4 + i] = DirichletMultinomialLogProbability(alphas, data);
	}
	double scale = result.maxCoeff();
	return (result - scale).exp();
}

HaploidProbs SequenceProb::HaploidSequencing(ReadData data) {
	HaploidProbs result;
	double alphas_total = (1.0 - params.phi_haploid) / params.phi_haploid;
	for (int i : { 0, 1, 2, 3 }) {
		double alphas[4];
		for (int k : { 0, 1, 2, 3 }) {
			if (k == i)
				alphas[k] = (1.0 - params.error_prob) * alphas_total;
			else
				alphas[k] = params.error_prob / 3.0 * alphas_total;
		}
		result[i] = DirichletMultinomialLogProbability(alphas, data);
	}
	double scale = result.maxCoeff();
	return (result - scale).exp();
}


void SequenceProb::UpdateLikelihood() {

//	auto it = site_data.all_reads.begin();
//	DiploidProbs anc_genotypes = DiploidSequencing(params, site_data.reference, *it);
//	anc_genotypes *= pop_genotypes;
	DiploidProbs anc_genotypes2 = anc_genotypes;
	DiploidProbs num_genotypes = anc_genotypes;
	for (auto t : allHap) {

		anc_genotypes2 *= (mutation.matrix() * t.matrix()).array();
		num_genotypes *= (non_mutation.matrix() * t.matrix()).array();

	}

	likelihood = 1- num_genotypes.sum() / anc_genotypes2.sum();
	//cerr << "\n" << anc_genotypes/anc_genotypes.sum() << endl;

	likelihood = anc_genotypes2.sum();
}

double SequenceProb::GetLikelihood() {
	return likelihood;
}

//ModelInput GetData();
ModelInput SequenceProb::GetData(){
	return data;
}
