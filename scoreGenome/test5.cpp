#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <sys/time.h>
#include <math.h>
//#include "CModel.h"
//#include "CSequences.h"
#include "scoreGenome.h"
using namespace std;

int main(int argc, char* argv[]) {
  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  
  cout << "USAGE: ./test5 [<modelFile> <genomeFile> <outFileName> <kMM-order> <substringSize>]" << endl;

  char * modelFileName = argv[1];
  char * genomeFileName = argv[2];
  char * outputFile = argv[3];
  int order = atoi(argv[4]);
  int substrLen = atoi(argv[5]);
  
  scoreGenome g;
  g.scoreGen(modelFileName, genomeFileName, outputFile, order, substrLen);
  cout << "Done" << endl;
  // m.printModel(outFileName);

	/*
  CSequences gen(genomeFileName);
  char * genome = gen.getGenome();

  char * modelFileName = argv[1];
  char * outFileName = argv[2];
  int order = atoi(argv[3]);
  
  CModel m(modelFileName, order);
  cout << "Finished reading in model" << endl;
  CModel m2(m); // create m2 of order 1 less than that of m 
  m2.printModel(outFileName);
  // m2.printModel(outFileName + ".5.param.txt");
  CModel m3(m2); // create m3 of order 1 less than that of m2
  // m3.printModel(outFileName + ".4.param.txt");
  CModel m4(m3); // create m3 of order 1 less than that of m2
  // m4.printModel(outFileName + ".3.param.txt");
  CModel m5(m4); // create m3 of order 1 less than that of m2
  // m5.printModel(outFileName + ".2.param.txt");
  strcpy(outFileName, "tmp.2.txt");
  m5.printModel(outFileName);
  */
  
  /*
  CModel m(7);
  
  char * str;
  str = new char[8];
  int len;
  strcpy(str, "AAACAGATAA");
  len = strlen(str);
  */
  /*
  for (len = 0; len < 15; len++) {
    str[len] = getchar();
    // cout << str[len] << endl;
    if (str[len] == '\n') {
      // cout << "hi" << endl;
      break;
    }
  }
  */
  /*
  cout << "String of length " << len << " read in" << endl;

  int pos;
  if (len >= 8) {pos = 8;} 
  else {pos = len;}
  int tmpMap = m.mapKmer(str,0,pos);
  cout << "Mapped index = " << tmpMap << endl;

  tmpMap = m.updateMapKmer(str,pos,pos,tmpMap);
  cout << "Next Mapped index - " << tmpMap << endl;
  tmpMap = m.updateMapKmer(str,pos+1,pos,tmpMap);
  cout << "Another Mapped index - " << tmpMap << endl;
  */
   
  // Now read the genome
  /* 
  string genomeFileName = "genomes/PA01.fasta";
  ifstream genomeFile;
  genomeFile.open(genomeFileName.c_str());
  if (!genomeFile.is_open()) {
    cout << "Could not open file " << genomeFileName << endl;
    return 0; 
  } 
	  
  char ch;
  string genome;
  getline(genomeFile, genome); 
  genome = ""; // discard first line
  
  while (genomeFile.get(ch)) {
    if (ch != '\n') {
      genome.push_back(ch);
    }
  }

  int len = genome.length();
  char* str;
  str = new char[len];
  strcpy(str, genome.c_str());
  for (int i = 0; i < 100; i++) {
    cout << str[i];
  }
  cout << endl << "Genome of length " << len << " read in" << endl;
  */

  gettimeofday(&tv2, NULL);
  double tm = (double) (tv2.tv_usec - tv1.tv_usec)/1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
  cout << "Time taken in execution = " << tm << " seconds\n";
  
  return 0;
} 
