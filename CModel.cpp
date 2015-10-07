/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or giri@cs.fiu.edu
 */

#include "CModel.h"
#include "CSequences.h"

CModel::CModel(int order) {
  this->order = order;
  this->size = pow(4, order+1) + pow(4, order); 
  this->logProbabilities = new float[size];
  // Do I need to initialize all log(Probabilities) to 0?
}

CModel::CModel(char * fName, int order) {
  // Model is constructed by reading a model file fName with order
  this->order = order;
  this->size = pow(4, order+1) + pow(4, order); 
  this->logProbabilities = new float[size];
  
  int index;
  float value;
  ifstream f;
  f.open(fName); // Open model file
  if (f.is_open()) {
    try { // In case there's something in the model's folder that shouldn't be there
      cout << "Reading Model from: " << fName << ":\n"; 
      int i = 0;
      while(f >> index >> value) {
	// cout << "Model value1: " << value << "\n";
	this->logProbabilities[i] = value; // Store model values
	i++;
	// cout << "Model value: " << value << "\n";
      }
    } catch (...) {}
  }

  // DEBUG test
  // for (int i = 0; i < size; i = i+4) {
    /*
    cout << this->logProbabilities.at(i) << "\t" << this->logProbabilities.at(i+1) << "\t" 
	 << this->logProbabilities.at(i+2) << "\t" << this->logProbabilities.at(i+3) << endl;
    */
  // } 
  // cout << "Finished reading model" << endl;
}

CModel::CModel(const CModel & big) {
  // Current model is constructed by copying m and decrementing its order
  int bigo =  big.order;
  this->order = bigo - 1;
  // int smallo = bigo - 1;
  int sPart1 = pow(4, order + 1); // size of first part of smaller model
  this->size = pow(4, order + 1) + pow(4, order); // size of smaller model
  // cout << bigo << "\t" << order << "\t" << sPart1 << "\t" << size << endl;
  this->logProbabilities = new float[size];
  // this->logProbabilities.assign(size,0.0); //Initialize all log(Probabilities) to 0

  // calculate probabilities for current model
  // initialize index to first part and 2nd part of logprobabilites vector
  int part1 = 0; // index to part storing conditional probabilites (part 1)
  int part2 = sPart1; // index to part storing initial probabiliites (part 2)

  // cout << "Starting to construct smaller model" << endl;
  for (int i = 4 * sPart1; i < 4 * size; i = i + 4) {
    double tmp = 0;
    for (int j = 0; j < 4; j++) {
      tmp += exp((double) big.logProbabilities[i+j]);
    }
    // cout << tmp << endl;
    double tmpLog = log(tmp);
    // cout << ": " << tmpLog << endl;;
    this->logProbabilities[part2] = (float) tmpLog;
    // cout << ":== " << tmpLog << endl;;
    part2++;
    for (int j = 0; j < 4; j++) {
      this->logProbabilities[part1+j] = (float) (big.logProbabilities[i+j] - tmpLog);
    }
    part1 += 4; 
    // if (part1 % 100 == 0) {cout << part1 << " "; }
    // if (part2 > 5010) {cout << part1 << " " << part2 << endl; }
  }
  // cout << "Finished constructing smaller model" << endl;
}

CModel::~CModel() {
  this->order=0;
  // this->frequency.clear();
  // this->logProbabilities.clear();
  // this->kmers->~CKmers();
  delete[] this->logProbabilities;
}

int CModel::mapKmer(char * str, int loc, int len) {
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
  if (len == this->order) {
    mapVal += pow(4,len+1);
  }
  return mapVal;
}

int CModel::updateMapKmer(char * str, int loc, int len, int prevVal) {
  // Assumes that the substring of length len starting from position loc-1 has map value prevVal
  // Now it computes the mapped value for substring of same length, but starting from location loc
  int mapVal = prevVal;
  int chVal;
  switch(str[loc]) {
  case 'A':
    chVal = 0; break;
  case 'T':
    chVal = 1; break;
  case 'G':
    chVal = 2; break;
  case 'C':
    chVal = 3; break;
  default:
    chVal = -1; break;
  }
  mapVal = 4 * mapVal;
  mapVal += chVal;

  switch(str[loc-len]) {
  case 'A':
    chVal = 0; break;
  case 'T':
    chVal = pow(4,len); break;
  case 'G':
    chVal = 2*pow(4,len); break;
  case 'C':
    chVal = 3*pow(4,len); break;
  default:
    return -1;
  }
  mapVal = mapVal - chVal;

  return mapVal;
}

void CModel::calcFreqsProbs(char * genome) {
  cout << "Calculating frequencies... " << (unsigned) strlen(genome) << "\n";

  int kmerStart = 0;
  int kmerEnd = order - 1;
  int genLength = (unsigned) strlen(genome);
  int * frequency = new int[this->size]; // temp array of frequencies
  for (int i = 0; i < this->size; i++) { frequency[i] = 1; } // initialize
  int nextN = -1; // keeps track of location of next 'N' in genome
  char * pch = strchr(genome, 'N'); // locate next 'N' and set nextN
  if (pch == NULL) { 
    nextN = genLength;
  } else {
    nextN = (int) (pch - genome);
  }
  cout << "Next occurrence of N in genome at " << nextN << "\t" << kmerEnd << "\t" << genLength << endl;

  // cout << ".";    
  int mappedIndex = -1; // used to map kmer to index in logProbabilities
  // cout << ".";    
  while (kmerEnd < genLength) {
    // first the freq of k-mers
    // cout << ".";    

    while (nextN <= kmerEnd) { // finds first kmer in genome without N's
      // if (kmerEnd > 6264400) {cout << "1:\t" << kmerEnd << "\t" << nextN << endl; }
      kmerStart = nextN + 1; // move window past nextN
      kmerEnd = kmerStart + order -1;
      if (kmerEnd  >= genLength) { break; } 
      pch = strchr(genome+kmerStart, 'N'); // locate next 'N and set nextN
      if (pch == NULL) { 
	nextN = genLength;
      } else {
	nextN = (int) (pch - genome);
      }
    }
    
    if (kmerEnd  >= genLength) { break; } 
    // cout << ".";
    // if (kmerEnd % 10000 == 0) {cout << "."; }
    // if (kmerEnd % 100000 == 0) {cout << kmerEnd << endl; }
    // Now we have a kmer with no 'N' between kmerStart and kmerEnd
    mappedIndex = mapKmer(genome, kmerStart, order); // get the Index of the k-mer
    frequency[mappedIndex]++; // incr frequency of k-mer
    // if (kmerEnd > 6264400) {cout << "2\t" << mappedIndex << "\t" << frequency[mappedIndex] << endl; }
    
    // next freq of (k+1)-mers
    while (nextN <= kmerEnd + 1) { // finds first (k+1)-mer in genome without N's
      if (kmerEnd > 6264400) {cout << "3::\t" << kmerEnd << "\t" << nextN << endl; }
      kmerStart = nextN + 1;
      kmerEnd = kmerStart + order -1;
      if (kmerEnd  >= genLength) { break; } 
      pch = strchr(genome+kmerStart, 'N'); // locate next 'N and set nextN
      if (pch == NULL) { 
	nextN = (unsigned) strlen(genome);
      } else {
	nextN = (int) (pch - genome);
      }
    }
    if (kmerEnd  >= genLength) { break; } 
    // if (kmerEnd % 10000 == 0) {cout << "."; }
    // if (kmerEnd % 100000 == 0) {cout << kmerEnd << endl; }
    mappedIndex = mapKmer(genome, kmerStart, order+1);
    frequency[mappedIndex]++; // incr frequency of (k+1)-mer
    // if (kmerEnd > 6264400) {cout << kmerEnd << "\t" << mappedIndex << endl; }
    // if (kmerEnd > 6264400) {cout << "4\t" << kmerEnd << "\t" << nextN << endl; }
    
    
    kmerStart++;
    kmerEnd++;
  }

  // unclear if this is really needed any more. Need to check.
  /* 
  // We need to add the last (k+1)-mer
  try {
    mappedIndex = this->kmers->getKmerList().at(tmpKmer);
    this->frequency.at(mappedIndex)++;
    //cout << "Last one!! " << tmpKmer << ": " << this->frequency.at(mappedIndex) << "\n";
  } catch(...) {} //If there's an N in the k-mer, just skip it
  */

  cout << "outside while loop" << endl;
  // Calculating all Probabilites from frequencies of (k+1)-mers
  int sum = 0; int totalSum = 0; 
  int length = pow(4,order+1); // first 4^{order+1} frequencies
  cout << "Calculating probabilities...\n";
  for (int i=0; i<=length-4; i+=4) {
    sum = frequency[i] + frequency[i+1] + frequency[i+2] + frequency[i+3];
    for(int j=i; j<i+4; j++) {
      this->logProbabilities[j] = log((float)frequency[j] / sum);
      // cout << i << "," << j << " = "<< this->logProbabilities[j] << "\n";
    }
  }
  for(int i=length; i<this->size; i++) {
    totalSum += frequency[i];
  }
  cout << totalSum << " total Sum " << endl;
  for(int i=length; i<this->size; i++) {
    this->logProbabilities[i] = log((float)frequency[i] / totalSum);
    // cout << i << " = "<< this->logProbabilities[i] << "\n";
  }
  delete[] frequency;
}


void CModel::buildModel(char * genome, char * outputName) {
  this->calcFreqsProbs(genome); //1. Obtain frequencies & probs

  // this->printModel(outputName); 
  //2. Write log(probabilities) to file
  cout << "finsihed calc freqs and probs" << endl; 
  ofstream outputFile;
  outputFile.open(outputName);

  if (outputFile.is_open()) {
    for(int i=0; i<this->size; i++) {
      outputFile << i << "\t" << this->logProbabilities[i] << "\n";
    }
  }
}

void CModel::printModel(char * outFileName) {
  ofstream f;
  f.open(outFileName);
  int lim = 5 * pow(4,order);
  for (int i = 0; i < lim; i = i+4) {
    f << this->logProbabilities[i] << "\t" << this->logProbabilities[i+1] << "\t" 
      << this->logProbabilities[i+2] << "\t" << this->logProbabilities[i+3] << endl;
  }
  f.close();
}
