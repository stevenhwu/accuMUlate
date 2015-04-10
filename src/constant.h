//
// Created by steven on 4/3/15.
//
#pragma once
#ifndef _CONSTANT_H_
#define _CONSTANT_H_



const int ANCESTOR_COUNT = 10;
const int BASE_COUNT = 4;

const Eigen::IOFormat nice_row(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");

const double ERROR_THRESHOLD = 1e-14;//std::numeric_limits<double>::epsilon()*10;

#endif //_CONSTANT_H_
