/*
Authors: Jiangrui Yu
Copyright:
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

* This code provides improvements over exact relu function. Following improvements are made:
* - Relu and Trunc Fusion: no extra overhead for Trunc after relu
* - ReluTrunc and Extend fusion: bitwidth varies in different layers
* - Trade off for number of pieces
* - fuse carry onto the next level to eliminate the end overhead
* - Eliminate overhead fo
r Auto's bit
*/

#ifndef RELUTRUNC_H__
#define RELUTRUNC_H__

#include "BuildingBlocks/aux-protocols.h"
#include "Millionaire/equality.h"
#include "Millionaire/millionaire_with_equality.h"

class ReluTrunc2 {
public:
  sci::IOPack *iopack;
  sci::OTPack *otpack;
  sci::PRG128 *prg;
  AuxProtocols *aux;
  int party;

  // Constructor
  ReluTrunc2(int party, sci::IOPack *iopack, sci::OTPack *otpack);

  // Destructor
  ~ReluTrunc2();

  // Relu and then Truncate (right-shift) by shifting in the same ring (round towards -inf)
  void ReluTrunc2::relutrunc(int32_t num,
                             uint64_t *input,
                             uint64_t *output,
                             int32_t total_bits,
                             int32_t auto_bits,
                             int32_t trunc_bits,
                             int32_t res_bits,
                             int32_t final_bits);
  void ReluTrunc2::wrap_computation(int32_t num,
                        uint64_t *input, uint8_t *output, int32_t bitlength);
};

#endif // RELUTRUNC_H__
