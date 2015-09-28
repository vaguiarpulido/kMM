/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or vaguiarpulido@gmail.com
 */

#include<stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <math.h>
using namespace std;

#ifndef CKMERS_H_
#define CKMERS_H_

class CKmers {
public:
	CKmers(int order);
	virtual ~CKmers();
	map<string, int> getKmerList();
private:
	int order;
	map<string, int> kmerList;

	char next(char current);
	void successor(string* current, int* pos);
	bool generateKmerList();
};

#endif /* CKMER_H_ */
