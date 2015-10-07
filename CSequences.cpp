/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or giri@cs.fiu.edu
 */

#include "CSequences.h"

CSequences::CSequences(char * inputFile) {
  this->input.open(inputFile);
}

/*
bool CSequences::LoadSequences()
{
  if (!this->input.is_open()) return false;
  string sequence="";
  char c;
  bool first=true;
  
  while (this->input.get(c)) {
    c = toupper(c); //Putting everything to upper case
    switch(c) {
    case 'A': case 'C': case 'G': case 'T': case 'N':
      sequence.push_back(c);
      break;
    case '>':
      if (first) {
	first = false;
      } else { //End of sequence
	if (!sequence.empty()) {
	  this->sequences.push_back(sequence);
	  //cout << "Sequence: " << sequence << "\n";
	  sequence.clear();
	}
      }
      while(this->input.get(c)) { //Skip lines with names
	if (c=='\n') break;
      }
      break;
    default:
      break;
    }
  }
  if (!sequence.empty()) { //Save last line
    this->sequences.push_back(sequence);
    //cout << "Sequence: " << sequence << "\n";
  }
  //cout << "Sequences: " << this->sequences.size() << "\n";  
  return true;
}
*/

CSequences::~CSequences() {
  if (this->input.is_open()) {
    input.close();
  }
  // this->sequences.clear();
}


char * CSequences::getGenome() {
  if (!input.is_open()) {
    cout << "Could not open genome file " << endl;
    return 0; 
  } 
	  
  char ch;
  string genome;
  getline(input, genome); 
  genome = ""; // discard first line
  while (input.get(ch)) {
    if (ch != '\n') {
      genome.push_back(ch);
    }
  }

  int len = genome.length();
  genomeSeq = new char[len];
  strcpy(genomeSeq, genome.c_str());
  for (int i = 0; i < 100; i++) {
    cout << genomeSeq[i];
  }
  cout << endl << "Genome of length " << len << " read in" << endl;

  return this->genomeSeq;
}
