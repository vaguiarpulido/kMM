/*
 * CReads.cpp
 *
 *  Created on: 21 de sept. de 2015
 *      Author: vanessa
 */

#include "CReads.h"

CReads::CReads(string inputFile) {
	this->input.open(inputFile.c_str());
}

bool CReads::LoadReads()
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

	//cout << "Reads: " << this->sequences.size() << "\n";

	return true;
}


CReads::~CReads() {
	if (this->input.is_open()) {
		input.close();
	}
	this->sequences.clear();
}


vector<string> CReads::getSequences() {
	if (this->sequences.empty()) { //First time initialization
		this->LoadReads();
	}

	return this->sequences;
}

/*
int main(int argc, char *argv[]) {
	CReads *reads = new CReads("test.txt");
	reads->getSequences();
	reads->~CReads();
}
*/
