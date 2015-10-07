/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or giri@fiu.edu
 */

#ifndef CSEQUENCES_H_
#define CSEQUENCES_H_

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
using namespace std;

class CSequences {
public:
	CSequences(char * inputFile);
	virtual ~CSequences();
	char * getGenome();

private:
	ifstream input;
	// vector<string> sequences;
	char * genomeSeq;
	// bool LoadSequences();
};

#endif /* CSEQUENCES_H_ */
