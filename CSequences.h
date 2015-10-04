/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or vaguiarpulido@gmail.com
 */

#ifndef CSEQUENCES_H_
#define CSEQUENCES_H_

#include <iostream>
#include<stdio.h>
#include <fstream>
#include <vector>
using namespace std;

class CSequences {
public:
	CSequences(string inputFile);
	virtual ~CSequences();
	vector<string> getSequences();

private:
	ifstream input;
	vector<string> sequences;

	bool LoadSequences();
};

#endif /* CSEQUENCES_H_ */
