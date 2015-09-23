/*
 * CReads.cpp
 *
 *  Created on: 21 de sept. de 2015
 *      Author: vanessa
 */

#include "CReads.h"

CReads::CReads(const char* inputFile) {
	this->input.open(inputFile);

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
					cout << "Sequence: " << sequence << "\n";
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
		cout << "Sequence: " << sequence << "\n";
	}

	cout << "Reads: " << this->sequences.size() << "\n";

	return true;

}

CReads::~CReads() {
	if (this->input.is_open()) {
		input.close();
	}
}

/*
int main(int argc, char *argv[]) {
	CReads *reads = new CReads("test.txt");
	reads->LoadReads();
	reads->~CReads();
}
*/
