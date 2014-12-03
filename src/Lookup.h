#pragma once
#ifndef __Lookup_H_
#define __Lookup_H_

#include "string"

namespace LookupTable {
//    class Lookup{
//struct LookupTable {
//
//    extern const int ANCESTOR_COUNT;
//    extern const int BASE_COUNT;


    extern const int index_vector[4][10]; //Same as summary stat now, number of mismatches between D (o) and A(m,f)
    extern const double summary_stat_index_lookup[4][10];

    extern const double summary_stat_same_lookup_table[10][4];

/*
10 cat version
  AA AC AG AT CC CG CT GG GT TT

  descnt = A
      A  C  G  T
 *AA == !! !! !!
  AC =! =! !! !!
  AG =! !! =! !!
  AT =! !! !! =!
 *CC !! == !! !!
  CG !! =! =! !!
  CT !! =! !! =!
 *GG !! !! == !!
  GT !! !! =! =!
 *TT !! !! !! ==


*/
    extern const int index_converter_16_to_10_single[16];
    extern const int index_converter_4_4_to_10[4][4];

    extern const int index_converter_4_4_to_16[4][4];
    extern const int index_converter_10_to_16[10];

    extern const std::string genotype_lookup_10[10];



};
#endif //__Lookup_H_
