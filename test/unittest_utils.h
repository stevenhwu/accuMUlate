//
// Created by steven on 4/3/15.
//

#ifndef _UNITTEST_UTILS_H_
#define _UNITTEST_UTILS_H_

#include "gtest/gtest.h"
#include "model.h"
#include "constant.h"


extern const double ERROR_THRESHOLD;

extern void ASSERT_GENOTYPES(HaploidProbs expected, HaploidProbs data);

extern void ASSERT_GENOTYPES(DiploidProbs expected, DiploidProbs data);


#endif //__UNITTEST_UTILS_H_
