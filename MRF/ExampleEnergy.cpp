#include "mrf.h"
#include "ICM.h"
#include "GCoptimization.h"
#include "MaxProdBP.h"
#include "TRW-S.h"
#include "BP-S.h"

#include <ctime>
#include <new>
#include <iostream>

const int sizeX = 50;
const int sizeY = 50;
const int numLabels = 20;

MRF::CostVal smoothnessCost(int node1, int node2, int label1, int label2);
MRF::CostVal dataCost(int node, int label);

int main(int argc, char **argv)
{

  double lowerBound;
  float t,tot_t;
  int iter;

  int seed = 1124285485;
  srand(seed);

  DataCost *data         = new DataCost(dataCost);
  SmoothnessCost *smoothness = new SmoothnessCost(smoothnessCost);
  EnergyFunction *energy = new EnergyFunction(data,smoothness);

  MRF* mrf = new MaxProdBP(sizeX,sizeY,numLabels,energy);
  mrf->initialize();
  mrf->clearAnswer();

  MRF::EnergyVal E = mrf->totalEnergy();
  std::cout << "Energy at the Start= " << (float)E << " (" << (float)mrf->smoothnessEnergy()
                                       << ", " << (float)mrf->dataEnergy() << ")" << std::endl;;

  tot_t = 0;
  for (iter=0; iter < 10; iter++)
  {
    mrf->optimize(1, t);

    E = mrf->totalEnergy();
    tot_t = tot_t + t ;
    std::cout << "energy = " << (float)E << " (" << tot_t << " secs)" << std::endl;
  }

  delete mrf;

  return 0;
}

MRF::CostVal dataCost(int pix, int i)
{
  return ((pix*i + i + pix) % 30) / ((MRF::CostVal) 3);
}

MRF::CostVal smoothnessCost(int node1, int node2, int label1, int label2)
{
  // The cost of assigning label1 to node1 and label2 to node2

  MRF::CostVal answer = (node1*(label1+1)*(label2+2) + node2*label1*label2*node1 - 2*label1*label2*node1) % 100;
  return answer / 10;
}
