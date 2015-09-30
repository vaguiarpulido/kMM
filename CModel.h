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
#include<stdio.h>
#include <vector>
#include <math.h>
#include "CKmers.h"
using namespace std;

class CModel {
public:
	CModel(int order);
	virtual ~CModel();
	void buildModel(string genome, string outputName);

private:
	int order;
	vector<int> frequency;
	vector<float> logProbabilities;
	CKmers* kmers;

	void calculateFrequencies(string genome);
	void calculateProbabilities();
};

#endif /* CMODEL_H_ */
