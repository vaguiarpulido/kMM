/*
 * CReads.h
 *
 *  Created on: 21 de sept. de 2015
 *      Author: vanessa
 */

#ifndef CREADS_H_
#define CREADS_H_

#include <iostream>
#include<stdio.h>
#include <fstream>
#include <vector>
using namespace std;

class CReads {
public:
	CReads(const char* inputFile);
	bool LoadReads();
	virtual ~CReads();
private:
	ifstream input;
	vector<string> sequences;
};

#endif /* CREADS_H_ */
