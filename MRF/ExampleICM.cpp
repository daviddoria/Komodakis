#include "mrf.h"
#include "ICM.h"
#include "GCoptimization.h"
#include "MaxProdBP.h"
#include "TRW-S.h"
#include "BP-S.h"
#include "ExampleEnergyFunctions.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>

const int sizeX = 50;
const int sizeY = 50;
const int numLabels = 20;

MRF::CostVal D[sizeX*sizeY*numLabels];
MRF::CostVal V[numLabels*numLabels];
MRF::CostVal hCue[sizeX*sizeY];
MRF::CostVal vCue[sizeX*sizeY];

#ifdef COUNT_TRUNCATIONS
int truncCnt, totalCnt;
#endif

MRF::CostVal fnCost(int pix1, int pix2, int i, int j);
MRF::CostVal dCost(int pix, int i);

int main(int argc, char **argv)
{
  MRF* mrf;
  EnergyFunction *energy;
  MRF::EnergyVal E;
  double lowerBound;
  float t,tot_t;
  int iter;

  int seed = 1124285485;
  srand(seed);

  int Etype = 0;

  energy = generate_DataARRAY_SmoothFIXED_FUNCTION();
  //energy = generate_DataARRAY_SmoothTRUNCATED_LINEAR();
  //energy = generate_DataARRAY_SmoothTRUNCATED_QUADRATIC();
  //energy = generate_DataFUNCTION_SmoothGENERAL_FUNCTION();


  printf("\n*******Started ICM *****\n");

  mrf = new ICM(sizeX,sizeY,numLabels,energy);
  mrf->initialize();
  mrf->clearAnswer();

  E = mrf->totalEnergy();
  printf("Energy at the Start= %g (%g,%g)\n", (float)E,
	  (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

  tot_t = 0;
  for (iter=0; iter<6; iter++) {
      mrf->optimize(10, t);

      E = mrf->totalEnergy();
      tot_t = tot_t + t ;
      printf("energy = %g (%f secs)\n", (float)E, tot_t);
  }

  delete mrf;


  return 0;
}

MRF::CostVal dCost(int pix, int i)
{
    return ((pix*i + i + pix) % 30) / ((MRF::CostVal) 3);
}

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{
    if (pix2 < pix1) { // ensure that fnCost(pix1, pix2, i, j) == fnCost(pix2, pix1, j, i)
	int tmp;
	tmp = pix1; pix1 = pix2; pix2 = tmp;
	tmp = i; i = j; j = tmp;
    }
    MRF::CostVal answer = (pix1*(i+1)*(j+2) + pix2*i*j*pix1 - 2*i*j*pix1) % 100;
    return answer / 10;
}