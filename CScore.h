/*
 * CScore.h
 *
 *  Created on: 23 de sept. de 2015
 *      Author: vanessa
 */

#include<stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "CReads.h"
#include "CKmers.h"

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
