/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or vaguiarpulido@gmail.com
 */

#include "CKmers.h"

CKmers::CKmers(int order) {
	this->order = order;
}

CKmers::~CKmers() {
	this->order = 0;
	this->kmerList.clear();
}

char CKmers::next(char current) {

	switch(current) {
	case 'A':
		return 'T';
	case 'T':
		return 'G';
	case 'G':
		return 'C';
	case 'C':
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
	int indexIniProb = pow(4, this->order+1);

	if (this->order<=0) return false;

	//Generate the first kmer (length k+1)
	for(int i=0; i<=kmerPos; i++) {
		current.push_back('A');
	}

	//Generate the rest of kmers
	while (kmerPos >= 0) {
		this->kmerList[current] = index;
		index++;
		if(index % 4 == 0) {
			this->kmerList[current.substr(0,this->order)] = indexIniProb;
			indexIniProb++;
		}
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
	CKmers *kmers = new CKmers(6);
	map<string,int> kmerMap = kmers->getKmerList();
	  // show content:
	  for (map<string,int>::iterator it=kmerMap.begin(); it!=kmerMap.end(); ++it)
	    std::cout << it->first << " \t " << it->second << '\n';
}
*/
