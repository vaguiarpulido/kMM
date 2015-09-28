/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or vaguiarpulido@gmail.com
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
	CReads(string inputFile);
	virtual ~CReads();
	vector<string> getSequences();

private:
	ifstream input;
	vector<string> sequences;

	bool LoadReads();
};

#endif /* CREADS_H_ */
