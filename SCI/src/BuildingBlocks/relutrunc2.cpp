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
#define OT_bits 7
#include "BuildingBlocks/relutrunc2.h"

using namespace std;
using namespace sci;

ReluTrunc2::ReluTrunc2(int party, IOPack *iopack, OTPack *otpack) {
  this->party = party;
  this->iopack = iopack;
  this->otpack = otpack;
  this->aux = new AuxProtocols(party, iopack, otpack);
  this->prg = new sci::PRG128;
}

ReluTrunc2::~ReluTrunc2() {
  delete this->aux;
  delete this->prg;
}

/*

 * The main function for relu.
 * The element will be split into the following parts:
 * |<--upper-bit-->|<--trunc-bit-->|<--auto-bit-->|
 * - auto bits are extra bits resulted from the automorphism operation
 *   - this part will be dropped directly
 * - trunc bits are the bit-width that will be truncated
 * - upper bits are the remaining bits
TODO: the function is not robust enough, could only handle legal cases
   - Invalid for 64 bits input
The protocol:
  1. compute the carry of the trunc bits, then compute the upper wrap and the msb of the upper bits
  Cost: 4 KKOT
  2. LUT the upper wrap with corresponding bit-width
  Cost: 1 IKNPOT
  3. LUT the lower wrap with corresponding bit-width
  Cost: 1 IKNPOT
  4. MUX
  Cost: 1IKNPOT
  total: 4KKOT+3IKNPOT ~440MB

*/
void ReluTrunc2::relutrunc(int32_t num, uint64_t *input, uint64_t *output, int32_t total_bits,
                          int32_t auto_bits, int32_t trunc_bits, int32_t upper_bits, int32_t final_bits) {

  uint64_t mask_all = (1ULL << (upper_bits+trunc_bits)) - 1;
  uint64_t mask_trunc = (1ULL << trunc_bits) - 1;
  uint64_t mask_upper = (1ULL << (upper_bits-1)) - 1; // minus one is for msb

  uint64_t *input_lower = new uint64_t[num];
  uint64_t *input_upper = new uint64_t[num];
  uint64_t *input_msb = new uint64_t[num];
  uint64_t *temp_output = new uint64_t[num];
  uint64_t *arith_wrap_total = new uint64_t[num];
  uint64_t *arith_wrap_lower = new uint64_t[num];
  uint8_t *wrap_lower = new uint8_t[num];
  uint8_t *wrap_upper = new uint8_t[num];
  uint8_t *drelu = new uint8_t[num];
  // We assume the last auto_bits is 0 as they are plaintext that donot need to be secret shared
  for (int i = 0; i < num;i++){
    input[i] = (input[i]>>auto_bits) & mask_all;
    input_lower[i] = (input[i] & mask_trunc);
    input_upper[i] = ((input[i]>>trunc_bits) & mask_upper); 
    input_msb[i] = (input[i]>>(upper_bits+trunc_bits-1)) & 1;
  }
  // We deal with the trunc and upper bits together. We need three bits: the carry of trunc part, the msb and upper-wrap of the upper part
  this->aux->wrap_computation(input_lower, wrap_lower, num, trunc_bits); 

  // We deal with the upper bits, the drelu result and the wrap_upper will be set into two bits
  if(party == sci::ALICE){

    uint64_t ** ot_messages=new uint64_t*[num];
    for (int i = 0; i < num;i++){
      //TODO: some random bits should be added
      ot_messages[i] = new uint64_t[(1<<upper_bits)];
      for(int j = 0; j < (1<<upper_bits); j++){
        // will wrap
        // j represnets: |msb|upper_share|lower_wrap(1 bit)|
        int upper_share = (j >> 1)&(mask_upper);
        int upper_msb_share = j>>(upper_bits);
        int lower_wrap_share = j & 1;
        // if lower bit will carry to the msb bit
        int carry_is_one = ((upper_share + input_upper[i]) > (1 << (upper_bits - 1))) || ((upper_share + input_upper[i]) == (1 << (upper_bits - 1))) && (static_cast<uint8_t>(lower_wrap_share) != wrap_lower[i]);
        // get the drelu share and wrap share
        int drelu = input_msb[i] ^ upper_msb_share ^ carry_is_one;
        int wrap = (input_msb[i] + upper_msb_share + carry_is_one >= 2);
        ot_messages[i][j] = (wrap<<1) + drelu;
      }
    }
    otpack->kkot[upper_bits]->send(ot_messages, num, 2);
  }else{
    uint8_t * ot_choices=new uint8_t[num];
    for (int i = 0;i< num;i++){
      ot_choices[i] = (input_upper[i]<<1)+wrap_lower[i];
    }
      otpack->kkot[upper_bits]->recv(wrap_upper, ot_choices, num, 2);
      for (int i = 0; i < num;i++){
        drelu[i] = wrap_upper[i] & 1;
        wrap_upper[i] = wrap_upper[i]>>1;
      }
  }
    this->aux->B2A(wrap_upper, arith_wrap_total, num, final_bits-upper_bits);
    this->aux->B2A(wrap_lower, arith_wrap_lower, num, final_bits);

  for (int i = 0; i < num; i++) {
    output[i]=((input[i]>>(auto_bits+trunc_bits)) + arith_wrap_lower[i] -
               (1ULL << (final_bits-upper_bits)) * arith_wrap_total[i]) & mask_all;
  }

    this->aux->multiplexer(drelu,temp_output,output,num,final_bits,final_bits);

  delete[] input_lower;
  delete[] input_upper;
  delete[] input_msb;
  delete[] temp_output;
  delete[] arith_wrap_total;
  delete[] arith_wrap_lower;
  delete[] wrap_lower;
  delete[] wrap_upper;
  delete[] drelu;
  return;
}

// A more comm-efficient wrap function
void ReluTrunc2::wrap_computation(int32_t num,
                      uint64_t *input, uint8_t *output, int32_t bitlength){

  int piece_length = OT_bits-1;
  int pieces = (bitlength+piece_length-1)/piece_length;
  int average_length = (bitlength+pieces-1)/pieces;
  uint64_t ** ot_messages=new uint64_t*[num];

  uint8_t* wrap =  new uint8_t[num]; 
  memset(output, 0, num);
  // encode the OT lut
  for (int i = 0; i < pieces;i++){
    int start = i*average_length;
    int end = (i+1)*average_length;
    if (end > bitlength){
      end = bitlength;
    }
    int length = end-start;
    uint64_t *input_piece = new uint64_t[num];
    for (int j = 0; j < num;j++){
      input_piece[j] = (input[j]>>start) & ((1ULL<<length)-1);
    }
    if(party == sci::ALICE){
    
    for (int j = 0; j < num;j++){
      ot_messages[j] = new uint64_t[(1<<(length+1))];
      for(int k = 0; k < (1<<(length+1)); k++){
        int upper_share = k >> 1;
        int lower_wrap_share = k & 1;
        int wrap_bit = (upper_share + input_piece[j] > (1 << (length))) || (upper_share + input_piece[j] == (1 << (length))) && (lower_wrap_share!=wrap[j]);
        ot_messages[j][k] = wrap_bit;
      }
    }
    otpack->kkot[length]->send(ot_messages, num, 1);
    }
  else{
    uint8_t * ot_choices=new uint8_t[num];
    for (int j = 0;j< num;j++){
      ot_choices[j] = input[j]<<1+wrap[j]; ;
    }
    otpack->kkot[length]->recv(wrap, ot_choices, num, 1);
  }
}
}