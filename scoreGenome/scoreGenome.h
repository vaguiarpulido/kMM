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
#include <string.h>
#include <vector>
#include <math.h>

using namespace std;

#ifndef scoreGenome_H_
#define scoreGenome_H_

class scoreGenome {
public:
	scoreGenome();
	virtual ~scoreGenome();
	void scoreGen(char * modelFileName, char * genomeFileName, char * outputFile, int order, int substrLen);

private:
	float scoreSubstr(char * genome, int loc, int length, int order, float * model);
	int mapKmer(char * str, int loc, int len, int order);
};

#endif /* scoreGenome_H_ */
