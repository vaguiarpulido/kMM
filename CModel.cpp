/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or vaguiarpulido@gmail.com
 */

#include "CModel.h"
#include <sys/time.h>
#include "CSequences.h"

CModel::CModel(int order) {

	this->order = order;
	this->kmers = new CKmers(order);
	int size = this->kmers->getKmerList().size();
	this->frequency.assign(size,1); //Initialize all frequencies to 1
	this->logProbabilities.assign(size,0.0); //Initialize all log(Probabilities) to 0
}

CModel::~CModel() {
	this->order=0;
	this->frequency.clear();
	this->logProbabilities.clear();
	this->kmers->~CKmers();
}

void CModel::calculateFrequencies(string genome) {
	cout << "Calculating frequencies... " << genome.length() << "\n";

	string tmpKmer = genome.substr(0, this->order); //From the beginning, take k characters
	int mappedIndex = -1;

	int i=this->order;
	while (i<genome.length()) {
		//k-mers
		try {
			mappedIndex = this->kmers->getKmerList().at(tmpKmer); //This will return the position of the k-mer
			this->frequency.at(mappedIndex)++; //We increase the frequency of the positions for initial probabilities
			//cout << tmpKmer << ": " << this->frequency.at(mappedIndex) << "\n";
		} catch (...) {} //If there's an N in the k-mer, just skip it

		//(k+1)-mers
		tmpKmer.push_back((char)genome.at(i));
		//cout << i << "\t" << tmpKmer;
		try {
			mappedIndex = this->kmers->getKmerList().at(tmpKmer);
			this->frequency.at(mappedIndex)++; //We increase the frequency of (k+1)-mers
			//cout << tmpKmer << ": " << this->frequency.at(mappedIndex) << "\n";
		} catch (...) {} //If there's an N in the (k+1)-mer, just skip it
		i++;
		tmpKmer.erase(tmpKmer.begin());
		//cout << "\t" << tmpKmer << "\n";
	}


	//We need to add the last (k+1)-mer
	try {
		mappedIndex = this->kmers->getKmerList().at(tmpKmer);
		this->frequency.at(mappedIndex)++;
		//cout << "Last one!! " << tmpKmer << ": " << this->frequency.at(mappedIndex) << "\n";
	} catch(...) {} //If there's an N in the k-mer, just skip it

}

void CModel::calculateProbabilities() {
	int sum = 0;
	int length = pow(4,this->order+1);
	cout << "Calculating probabilities...\n";
	  for (int i=0; i<=length-4; i+=4) {
		  sum = this->frequency.at(i) + this->frequency.at(i+1) + this->frequency.at(i+2) + this->frequency.at(i+3);
		  for(int j=i; j<i+4; j++) {
			  this->logProbabilities.at(j) = log((float)this->frequency.at(j) / sum);
			  //cout << i << "," << j << " = "<< this->logProbabilities.at(j) << "\n";
		  }
	  }
	  for(int i=length; i<this->frequency.size(); i++) {
		  this->logProbabilities.at(i) = log((float) this->frequency.at(i) / (pow(4,this->order) - this->order));
	  }

}


void CModel::buildModel(string genome, string outputName) {

	this->calculateFrequencies(genome); //1. Obtain the frequencies
	this->calculateProbabilities();	//2. Obtain the log(probabilities)

	//3. Write the log(probabilities) to a file
	ofstream outputFile;
	outputFile.open(outputName.c_str());

	if (outputFile.is_open()) {
		for(int i=0; i<this->logProbabilities.size(); i++) {
			outputFile << i << "\t" << this->logProbabilities.at(i) << "\n";
		}
	}
}


/*

int main(int argc, char* argv[]) {
	CModel* model = new CModel(6);

	struct timeval tv1, tv2;
	gettimeofday(&tv1, NULL);

	CSequences* genomes = new CSequences("NC_018265.fasta");
	//CSequences* genomes = new CSequences("test3.fa");
	cout << "Loading genome...\n";
	genomes->getSequences();

	cout << "Building model...\n";
	for(int i=0; i<genomes->getSequences().size(); i++) {
		cout << "Sequence: " << i << "\n";
		model->buildModel((string)genomes->getSequences().at(i),"testModel.txt");
	}

	gettimeofday(&tv2, NULL);
	double tm = (double) (tv2.tv_usec - tv1.tv_usec)/1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
	cout << "Time taken in execution = " << tm << " seconds\n";

	genomes->~CSequences();
	model->~CModel();



return 0;
}
*/

