//
// Created by steven on 4/3/15.
//

#include "unittest_utils.h"

void ASSERT_GENOTYPES(HaploidProbs expected, HaploidProbs data) {
    for (int i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], data[i], ERROR_THRESHOLD );
    }
}

void ASSERT_GENOTYPES(DiploidProbs expected, DiploidProbs data) {
    for (int i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], data[i], ERROR_THRESHOLD );
    }
}

