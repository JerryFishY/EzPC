/*
Authors: Jiangrui Yu
Copyright:
Copyright (c) 2021 Microsoft Research
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "BuildingBlocks/relutrunc.h"
#include "BuildingBlocks/value-extension.h"

using namespace std;
using namespace sci;

ReluTrunc::ReluTrunc(int party, IOPack *iopack, OTPack *otpack) {
  this->party = party;
  this->iopack = iopack;
  this->otpack = otpack;
  this->aux = new AuxProtocols(party, iopack, otpack);
  this->mill = this->aux->mill;
  this->mill_eq = new MillionaireWithEquality(party, iopack, otpack);
  this->eq = new Equality(party, iopack, otpack);
  this->triple_gen = this->mill->triple_gen;
}

ReluTrunc::~ReluTrunc() {
  delete this->aux;
  delete this->mill_eq;
  delete this->eq;
}


void ReluTrunc::relutrunc(int32_t dim, uint64_t *inA, uint64_t *outB,
                          int32_t shift, int32_t bw) {
  if (shift == 0) {
    memcpy(outB, inA, sizeof(uint64_t) * dim);
    return;
  }
  assert((bw - shift) > 0 && "Truncation shouldn't truncate the full bitwidth");
  assert(bw <= 64 && "Total bitwidth should be less than 64");
  assert(shift <= 64 && "Shift bitwidth should be less than 64");
  assert(inA != outB);

  uint64_t mask_all = (bw == 64 ? -1 : ((1ULL << bw) - 1));
  uint64_t mask_lower = (shift == 64 ? -1 : ((1ULL << shift) - 1));
  uint64_t mask_upper =
      ((bw - shift - 1) == 64 ? -1 : ((1ULL << (bw - shift - 1)) - 1));
  uint64_t mask_msbupper =
      ((bw - shift) == 64 ? -1 : ((1ULL << (bw - shift)) - 1));    

  // Since it's relutrunc, logical shift is enough
  uint64_t *inA_lower = new uint64_t[dim];
  uint64_t *inA_upper = new uint64_t[dim];
  uint64_t *inA_msbupper = new uint64_t[dim];
  uint64_t *outA = new uint64_t[dim];

  uint8_t *wrap_lower = new uint8_t[dim];
  uint8_t *wrap_upper = new uint8_t[dim];
  uint8_t *eq_upper = new uint8_t[dim];
  uint8_t *wrap_total = new uint8_t[dim];
  uint8_t *and_upper = new uint8_t[dim];
  uint8_t *drelu = new uint8_t[dim];
  for (int i = 0; i < dim; i++)
  {
      inA_lower[i] = inA[i] & mask_lower;
      inA_upper[i] = (inA[i] >> shift) & mask_upper;
      inA_msbupper[i] = (inA[i] >> shift) & mask_msbupper;
      if (party == BOB)
      {
          inA_upper[i] = (mask_upper - inA_upper[i]) & mask_upper;
      }
  }

  this->aux->wrap_computation(inA_lower, wrap_lower, dim, shift); // Wrap Lower

  this->mill_eq->compare_with_eq(wrap_upper, eq_upper, inA_upper, dim, (bw - shift - 1)); // Wrap upper

  this->aux->AND(wrap_lower, eq_upper, and_upper, dim);
  for (int i = 0; i < dim; i++)
    {
        wrap_upper[i] ^= and_upper[i];
        if(party == BOB){
            drelu[i] = (inA[i] >> (bw - 1)) ^ wrap_upper[i]^1;
        }else{
            drelu[i] = (inA[i] >> (bw - 1)) ^ wrap_upper[i];
        }
    }

    if (party == sci::ALICE)
    {
      PRG128 prg;
      prg.random_bool((bool *)wrap_total, dim);
      uint8_t **spec = new uint8_t *[dim];
      for (int i = 0; i < dim; i++)
      {
        spec[i] = new uint8_t[2];
        spec[i][0] = (inA[i] >> (bw - 1)) ^ wrap_total[i];
        spec[i][1] = 1 ^ wrap_total[i];
      }
      this->aux->lookup_table<uint8_t>(spec, nullptr, nullptr, dim, 1, 1);

      for (int i = 0; i < dim; i++)
        delete[] spec[i];
      delete[] spec;
    }
    else
    { // party == sci::BOB
      uint8_t *lut_in = new uint8_t[dim];
      for (int i = 0; i < dim; i++)
      {
        lut_in[i] = (inA[i] >> (bw - 1));
      }
      this->aux->lookup_table<uint8_t>(nullptr, lut_in, wrap_total, dim, 1, 1);

      delete[] lut_in;
    }

    uint64_t *arith_wrap_total = new uint64_t[dim];
    uint64_t *arith_wrap_lower = new uint64_t[dim];
    this->aux->B2A(wrap_total, arith_wrap_total, dim, shift);
    this->aux->B2A(wrap_lower, arith_wrap_lower, dim, bw);

  for (int i = 0; i < dim; i++) {
    outA[i]=(inA_msbupper[i] + arith_wrap_lower[i] -
               (1ULL << (bw - shift)) * arith_wrap_total[i]) & mask_all;
    }

    this->aux->multiplexer(drelu,outA,outB,dim,bw,bw);

  delete[] outA;
  delete[] inA_lower;
  delete[] inA_upper;
  delete[] wrap_lower;
  delete[] wrap_upper;
  delete[] eq_upper;
  delete[] and_upper;
  delete[] arith_wrap_lower;
  delete[] arith_wrap_total;

  return;
}
