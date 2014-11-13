/*
 * MutationProb.cpp
 *
 *  Created on: Nov 8, 2014
 *      Author: Steven Wu
 */

#include <MutationProb.h>
#include <cmath>
#include <iostream>

MutationProb::MutationProb(const ModelParams &model_params) {

	params = model_params; //TODO check copy method/ref

	beta0 = 1.0;
	for (auto d : params.nuc_freq) {
		beta0 -= d * d;
	}
	beta0 = 1.0 / beta0;
	UpdateMu();
}

MutationProb::~MutationProb() {
	// TODO Auto-generated destructor stub
}


void MutationProb::UpdateMu() {
	UpdateMu(params.mutation_rate);
}

void MutationProb::UpdateMu(double mu) {
	params.mutation_rate = mu;

	mutation = MutationAccumulation(params, false);
	MutationMatrix mt = MutationAccumulation(params, true);

	non_mutation = mutation - mt;

//	std::cout << mutation << endl;
//	std::cout << non_mutation << endl;
//	mutation -= 1;
//	mt = MutationMatrix::Zero();
//	std::cout << mutation << endl;
//	std::cout << non_mutation << endl;

	MutationAccumulation2(false);
//	UpdateLikelihood();
}

MutationMatrix MutationProb::MutationAccumulation2(bool and_mut) {

	using namespace Eigen;
	double beta = CalculateBeta();
//	printf("beta:%.10f\n", beta); //0.9999999852
	Eigen::Matrix4d m;
	Eigen::Matrix4d m2 = Eigen::Matrix4d::Zero();


	Matrix3d m3 = Matrix3d::Identity();
	m3.row(1) = Vector3d(4,5,6);
//	std::cout << m3 << std::endl;

	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			m(i, j) = params.nuc_freq[i] * (1.0 - beta);
		}
		m(i, i) += beta;
	}
//	element wise operation
	for (int i : { 0, 1, 2, 3 }) {
		double c = params.nuc_freq[i] * (1.0 - beta);
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
//	printf("beta:%.10f\n", beta); //0.9999999852
	Eigen::Matrix4d m;
	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			m(i, j) = params.nuc_freq[i] * (1.0 - beta);
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

	double beta = exp(-beta0 * params.mutation_rate);
	return beta;
}
