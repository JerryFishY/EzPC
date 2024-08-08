/*
Author(s): Jiangrui Yu
*/

#include "BuildingBlocks/relutrunc.h"
#include <iostream>

using namespace sci;
using namespace std;

int dim = 1024;
int bitlength = 12;
int sf = 8;
uint64_t mask = (bitlength == 64 ? -1 : ((1ULL << bitlength) - 1));

// vars
int party, port = 32000;
IOPack *iopack;
OTPack *otpack;
ReluTrunc *relutrunc;
PRG128 prg;

void fusion() {
  uint64_t *inA = new uint64_t[dim];
  uint64_t *outB = new uint64_t[dim];

  prg.random_data(inA, dim * sizeof(uint64_t));
  prg.random_data(outB, dim * sizeof(uint64_t));

  for (int i = 0; i < dim; i++) {
    inA[i] &= mask;
    outB[i] = 0;
  }

  uint8_t *msbA = nullptr;

  uint64_t num_rounds = iopack->get_rounds();
  relutrunc->relutrunc(dim, inA, outB, sf,bitlength);
  num_rounds = iopack->get_rounds() - num_rounds;
  cout << "Num rounds (Relu Truncation Fusion): " << num_rounds << endl;

  if (party == ALICE) {
    uint64_t *inA_bob = new uint64_t[dim];
    uint64_t *outB_bob = new uint64_t[dim];
    iopack->io->recv_data(inA_bob, sizeof(uint64_t) * dim);
    iopack->io->recv_data(outB_bob, sizeof(uint64_t) * dim);
    for (int i = 0; i < dim; i++) {
      inA[i] = (inA[i] + inA_bob[i]) & mask;
      outB[i] = (outB[i] + outB_bob[i]) & mask;
    }
    cout << "Testing for correctness..." << endl;
    for (int i = 0; i < dim; i++) {
      // cout << inA[i] << " " << outB[i] << endl;
      if(inA[i]>=(1<<(bitlength-1))){
        inA[i] = 0;
      }else{
        inA[i] = inA[i] >> sf;
      }
      assert(inA[i] == outB[i]);
    }
    cout << "Correct!" << endl;
  } else { // BOB
    iopack->io->send_data(inA, sizeof(uint64_t) * dim);
    iopack->io->send_data(outB, sizeof(uint64_t) * dim);
  }
}

int main(int argc, char **argv) {
  ArgMapping amap;
  amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
  amap.arg("p", port, "Port Number");
  amap.arg("d", dim, "Size of vector");
  amap.parse(argc, argv);

  iopack = new IOPack(party, port, "127.0.0.1");
  otpack = new OTPack(iopack, party);
  relutrunc = new ReluTrunc(party, iopack, otpack);

  cout << "<><><><> Relu And Truncation Fusion <><><><>" << endl;
  fusion();
}
