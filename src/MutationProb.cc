/*
 * MutationProb.cpp
 *
 *  Created on: Nov 8, 2014
 *      Author: Steven Wu
 */

#include "MutationProb.h"
#include <cmath>
#include <iostream>

#include "SequenceProb.h"

MutationProb::MutationProb(const ModelParams &model_params) {

	params = model_params; //TODO check copy method/ref

	beta0 = 1.0;
	for (auto d : model_params.nuc_freq) {
		beta0 -= d * d;
	}
	beta0 = 1.0 / beta0;

    for (int i = 0; i < 4; ++i) {
        frequency_prior[i] = model_params.nuc_freq[i];
//        beta0 += params.nuc_freq[i] * params.nuc_freq[i];
        for (int j = i; j < 4; ++j) {
            int index10 = index_converter_16_to_10[i][j];
            ancestor_prior[index10] = model_params.nuc_freq[i] * model_params.nuc_freq[j];
            if(i != j){
                ancestor_prior[index10] *= 2; //Count both AC and CA
            }
        }
    }

    UpdateMu(model_params.mutation_rate);
}

MutationProb::~MutationProb() {
	// TODO Auto-generated destructor stub
}


//void MutationProb::UpdateMu() {
//
//}

void MutationProb::UpdateMu(double mu0) {
	mutation_rate.mu = mu0;
    mutation_rate.one_minus_mu = 1 - mutation_rate.mu;
    CalculateBeta();
	mutation = MutationAccumulation(params, false);
	MutationMatrix mt = MutationAccumulation(params, true);
//
//	non_mutation = mutation - mt;

//	std::cout << mutation << endl;
//	std::cout << non_mutation << endl;
//	mutation -= 1;
//	mt = MutationMatrix::Zero();
//	std::cout << mutation << endl;
//	std::cout << non_mutation << endl;

	MutationAccumulation2(false);//TODO: change to this version later
//	UpdateLikelihood();
}

MutationMatrix MutationProb::MutationAccumulation2(bool and_mut) {

	using namespace Eigen;
	double beta = CalculateBeta();
//	printf("beta0:%.10f\n", beta0); //0.9999999852
	Eigen::Matrix4d m;
	Eigen::Matrix4d m2 = Eigen::Matrix4d::Zero();


	Matrix3d m3 = Matrix3d::Identity();
	m3.row(1) = Vector3d(4,5,6);
//	std::cout << m3 << std::endl;

	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			m(i, j) = frequency_prior[i] * (1.0 - beta);
		}
		m(i, i) += beta;
	}
//	element wise operation
	for (int i : { 0, 1, 2, 3 }) {
		double c = frequency_prior[i] * (1.0 - beta);
//		c = (i+1) * 1.1;
		m2.row(i) = Vector4d::Constant(c);
	}
	m2.diagonal() += Vector4d::Constant(beta);
//	std::cout << m << std::endl;
//	std::cout << "\n" << m2 << std::endl;

	//cerr << m << endl;
	MutationMatrix result;
	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			for (int k : { 0, 1, 2, 3 }) {
				result(i * 4 + j, k) = 0.0;
				if (!and_mut || i != k)
					result(i * 4 + j, k) += 0.5 * m(i, k);
				if (!and_mut || j != k)
					result(i * 4 + j, k) += 0.5 * m(j, k);
			}
		}
	}

	MutationMatrix result2 = MutationMatrix::Zero();;
		for (int i : { 0, 1, 2, 3 }) {
			for (int j : { 0, 1, 2, 3 }) {
				for (int k : { 0, 1, 2, 3 }) {
//					result(i * 4 + j, k) = 0.0;
					if (!and_mut || i != k)
						result2(i * 4 + j, k) += 0.5 * m2(i, k);
					if (!and_mut || j != k)
						result2(i * 4 + j, k) += 0.5 * m2(j, k);
				}
			}
		}
//	printf("%d\n", and_mut);
//	std::cout << result << std::endl;
//	std::cout << "\n" << std::endl;

//	std::cout << (result2 - result).sum() << "<-Should be 0" << std::endl;



	return result;
}

MutationMatrix MutationProb::MutationAccumulation(const ModelParams &params, bool and_mut) {

	double beta = CalculateBeta();
//	printf("beta0:%.10f\n", beta0); //0.9999999852
	Eigen::Matrix4d m;
	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			m(i, j) = frequency_prior[i] * (1.0 - beta);
		}
		m(i, i) += beta;
	}
//	std::cout << m << std::endl;

	//cerr << m << endl;
	MutationMatrix result;
	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			for (int k : { 0, 1, 2, 3 }) {
				result(i * 4 + j, k) = 0.0;
				if (!and_mut || i != k)
					result(i * 4 + j, k) += 0.5 * m(i, k);
				if (!and_mut || j != k)
					result(i * 4 + j, k) += 0.5 * m(j, k);
			}
		}
	}
//	printf("%d\n", and_mut);
//	std::cout << result << std::endl;

	return result;
}

MutationMatrix MutationProb::GetMutation() {

	return mutation;
}

MutationMatrix MutationProb::GetNonMutation() {
	return non_mutation;
}

double MutationProb::CalculateBeta(){

	beta = exp(-beta0 * mutation_rate.mu);
	return beta;
}

MutationRate MutationProb::GetMutationRate() {
    return mutation_rate;
}

double MutationProb::GetBeta() {
    return beta;
}

Array4D MutationProb::GetFrequencyPrior() {
    return frequency_prior;
}

Array10D MutationProb::GetAncestorPrior() {
    return ancestor_prior;
}
