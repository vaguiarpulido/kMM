/*
 *      Author: 	Giri Narasimhan
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or giri@fiu.edu
 */

#include "scoreGenome.h"
#include <sys/time.h>

scoreGenome::scoreGenome() {
}

scoreGenome::~scoreGenome() {
}

int scoreGenome::mapKmer(char * str, int loc, int len, int order) {
// Takes a substring of str starting from loc of length len and computes its mapped value
  int mapVal = 0;
  int chVal;
  
  for(int i = loc; i < loc+len; i++) {
    switch(str[i]) {
    case 'A':
      chVal = 0; break;
    case 'T':
      chVal = 1; break;
    case 'G':
      chVal = 2; break;
    case 'C':
      chVal = 3; break;
    default:
      return -1;
    }
    mapVal = 4 * mapVal;
    mapVal += chVal;
    // cout << "Mapped index = " << chVal << " " << mapVal << endl;
  }
  if (len == order) {
    mapVal += pow(4,len+1);
  }
  return mapVal;
}

float scoreGenome::scoreSubstr(char * genome, int loc, int length, int order, float * model) {
  char * pch = strchr(genome+loc, 'N');
  if (pch != NULL) { 
    if (pch - (genome + loc) <= length) {
      return 0.0;
    }
  }
  float score;
  int mappedIndex = -1;
  
  //We need to calculate the initial probabilities
  mappedIndex = mapKmer(genome, loc, order, order);
  score = model[mappedIndex];
  
  for(int j = loc; j < loc + length - order; j++) { //Calculate the score of a substring
    mappedIndex = mapKmer(genome, j, order+1, order); //This will return the position of the kmer
    score += model[mappedIndex];
  }
  //cout << "\nScore for sequence: " << score << "\n";
  
  return score;
}

void scoreGenome::scoreGen(char * modelFileName, char * genomeFileName, char * outputFile, int order, int substrLen) {
  
  struct timeval tvA, tvB, tvC;

  ifstream modelFile, genomeFile;
  ofstream scoreResults;
  int size = pow(4,order+1) + pow(4,order);
  float * model = new float[size];
  float value = 0.0;
  int index=-1, mappedIndex=-1;
  float tmpScore=0.0;
  float * scores;
  int len = 0;

  modelFile.open(modelFileName); //Open the file that contains the probabilities
  if (modelFile.is_open()) {
    cout << "Model: " << modelFileName << "\n";
    while(modelFile >> index >> value) {
      model[index] = value; //Store the model values
      //cout << "Model value: " << value << "\n";
    }
    
    // Now read the genome
    genomeFile.open(genomeFileName);
    if (!genomeFile.is_open()) {
      cout << "Could not open file " << genomeFileName << endl;
      return; 
    } 
    cout << "Finsihed reading model " << modelFileName << endl;

    char ch;
    string str;
    getline(genomeFile, str); // skip first line
    str = "";
    while (genomeFile.get(ch)) {
      if (ch != '\n') {
	str.push_back(ch);
      }
    }
    cout << "Length of " << modelFileName << " = " << str.length() << "\n";
	    
    len = str.length();
    char * genomeSeq = new char[len];
    strcpy(genomeSeq, str.c_str());
    for (int i = 0; i < 100; i++) {
      cout << genomeSeq[i];
    }
    cout << endl << "Genome of length " << len << " read in" << endl;
	
    scores = new float[len - substrLen + 1];
    for (int i = 0; i <= len - substrLen; i++) {
      scores[i] = this->scoreSubstr(genomeSeq, i, substrLen, order, model);
      // cout << "DEBUG: " << i << endl;
    }
    
    modelFile.close();
  } //End if model was loaded
    
  //Write the final scores to a file
  try {
    scoreResults.open(outputFile);
    // cout << "DEBUG: ready to output score results to " << this->scores.size() << outputFile << endl;
  } catch (...) {
    cout << "DEBUG: Problem opening file " << outputFile << endl;
  }
  if (scoreResults.is_open()) {
    // scoreResults << "Best score\tBest model\n";
    for (int i=0; i <= len - substrLen; i++) {
      scoreResults << i << "\t" << scores[i] << endl;
      // cout << "DEBUG: " << i << endl;
    }
    scoreResults.close();
  } else {
    cout << "Could not open file " << outputFile << endl;
  }
}
