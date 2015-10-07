/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or vaguiarpulido@gmail.com
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "CKmers.h"
#include "CSequences.h"

using namespace std;

#ifndef CSCORE_H_
#define CSCORE_H_

class CScore {
public:
	CScore();
	virtual ~CScore();
	void scoreModels(string modelsPath, string readsFileName, string outputFile, int order);

private:
	vector<float> scores;
	vector<string> modelNames;
	map<string, int> kmerMap;

	float scoreRead(string read, int order, map<string, int> kmerMap, vector<float> model);
};

#endif /* CSCORE_H_ */
