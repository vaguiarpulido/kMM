#include "CKmers.h"

CKmers::CKmers(int size) {
	this->order = size;
}

CKmers::~CKmers() {
	this->order = 0;
	this->kmerList.clear();
}

char CKmers::next(char current) {

	switch(current) {
	case 'A':
		return 'C';
	case 'C':
		return 'G';
	case 'G':
		return 'T';
	case 'T':
		return 'A';
	default:
		return 'N';
	}

}

void CKmers::successor(string* current, int* pos) {
	char next;

	next = this->next((char)(*current)[*pos]);

	while (next == 'A') {
		(*current)[*pos] = next;
		*pos = *pos-1;

		//If we reach the first position and the next is 'A' then we're done
		if(*pos==-1) break;

		next = this->next((char)(*current)[*pos]);
	}

	if (*pos >= 0) { //To control the case in which we're done
		(*current)[*pos] = next;
		*pos=this->order;
	}

}

bool CKmers::generateKmerList() {
	int kmerPos = this->order; //The kmers will have length k+1
	string current = "";
	int index = 0;

	if (this->order<=0) return false;

	//Generate the first kmer (length k+1)
	for(int i=0; i<=kmerPos; i++) {
		current.push_back('A');
	}

	//Generate the rest of kmers
	while (kmerPos >= 0) {
		this->kmerList[current] = index;
		index++;
		//cout << "Element " << index << ":" << current << "\n";
		this->successor(&current,&kmerPos);
	}

	return true;
}

map<string, int> CKmers::getKmerList() {
	//Initialize it if it hasn't been already done
	if (this->kmerList.empty()) {
		this->generateKmerList();
	}
	return this->kmerList;
}

/*
int main(int argc, char *argv[]) {
	CKmers *kmers = new CKmers(1);
	//kmers->generateKmerList();
}
*/
