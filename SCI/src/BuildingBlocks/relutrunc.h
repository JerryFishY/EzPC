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

* This code provides relu and trunc fusion, allowing for an exact truncation without extra overhead. 
*/

#ifndef RELUTRUNC_H__
#define RELUTRUNC_H__

#include "BuildingBlocks/aux-protocols.h"
#include "Millionaire/equality.h"
#include "Millionaire/millionaire_with_equality.h"

class ReluTrunc {
public:
  sci::IOPack *iopack;
  sci::OTPack *otpack;
  TripleGenerator *triple_gen;
  MillionaireProtocol *mill;
  MillionaireWithEquality *mill_eq;
  Equality *eq;
  AuxProtocols *aux;
  int party;

  // Constructor
  ReluTrunc(int party, sci::IOPack *iopack, sci::OTPack *otpack);

  // Destructor
  ~ReluTrunc();

  // Relu and then Truncate (right-shift) by shifting in the same ring (round towards -inf)
  void relutrunc(
      // Size of vector
      int32_t dim,
      // input vector
      uint64_t *inA,
      // output vector
      uint64_t *outB,
      // right shift amount
      int32_t shift,
      // Input and output bitwidth
      int32_t bw      
      );

};

#endif // RELUTRUNC_H__
