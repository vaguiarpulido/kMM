/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or vaguiarpulido@gmail.com
 */

#ifndef CMODEL_H_
#define CMODEL_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <sys/time.h>
#include "CKmers.h"
using namespace std; 

class CModel {
public:
  CModel(int order);
  CModel(char * fName, int order); // read model from file fName
  CModel(const CModel & m); // copy model m and decrement order by 1
  virtual ~CModel();

  int size; // size of logProbabilities = 4^{k+1} + 4^k

  void buildModel(char * genome, char * outputName);
  int mapKmer(char * str, int loc, int len);
  int updateMapKmer(char * str, int loc, int len, int prevVal);
  void printModel(char * outFileName);

private:
  int order;
  float * logProbabilities;

  void calcFreqsProbs(char * genome);
};

#endif /* CMODEL_H_ */
