/*
 *      Author: 	Vanessa Aguiar-Pulido
 *      			Postdoctoral Research Associate
 *      			Bioinformatics Research Group (BioRG)
 *      			Florida International University (FIU)
 *      			Miami, FL, USA
 *
 *     	Contact: 	vaguiarp@fiu.edu or vaguiarpulido@gmail.com
 */

#include "CScore.h"
#include <sys/time.h>

CScore::CScore() {


}

CScore::~CScore() {

}

float CScore::scoreRead(string read, int order, map<string, int> kmerMap, vector<float> model) {
	float score = 0.0;
	bool first = true;
	int mappedIndex = -1;
	string tmpKmer = read.substr(0, order); //From the beginning, take k characters;

	//We need to calculate the initial probabilities
	while(first) {
		try {
			mappedIndex = kmerMap.at(tmpKmer);
			score += model.at(mappedIndex);
			first = false;
		} catch (...) { //If there is an N in the initial k-mer, iterate past it
			tmpKmer.erase(tmpKmer.begin());
			tmpKmer.push_back((char)read.at(order));
			order++;
		}
	}

	//cout << "tmpKmer " << order << ": " << tmpKmer << "\n";

	tmpKmer.push_back((char)read.at(order)); //First (k+1)-mer

	for(int j=order+1; j<read.length(); j++) { //Calculate the score of a read
		//cout << "tmpKmer " << j << ": " << tmpKmer << "\n";
		try {
			mappedIndex = kmerMap.at(tmpKmer); //This will return the position of the kmer
			score += model.at(mappedIndex);
			//cout << "Partial score "<< j << ": " << score << "\n";
		} catch (...) {} //If there's an N, just skip it
		tmpKmer.erase(tmpKmer.begin());
		tmpKmer.push_back((char)read.at(j));
	}

	//cout << "tmpKmer " << read.length() << ": " << tmpKmer << "\n";

	//We need to add the last kmer
	try {
		mappedIndex = kmerMap.at(tmpKmer);
		score += model.at(mappedIndex);
	} catch (...) {} //If there's an N, just skip it

	//cout << "Score for read "<< read << ": " << score << "\n";

	return score;

}

void CScore::scoreModels(string modelsPath, string readsFileName, string outputFile, int order) {
	ifstream listFile, modelFile;
	ofstream scoreResults;
	string modelName="", modelFull="";
	vector<float> model;
	float value = 0.0;
	int index=-1, mappedIndex=-1;
	string tmpKmer="", tmpRead="";
	float tmpScore=0.0;
	//Prepare to load the reads
	CReads* reads = new CReads(readsFileName);
	//Prepare to get the list of possible kmers for a model
	CKmers* kmers = new CKmers(order);


	//Get the full list of models
	string command = "ls "+ modelsPath +" > models.txt";
	system(command.c_str());

	//Open the file containing the names of the models
	listFile.open("models.txt");
	//cout << "Let's open the models file\n";

	if (listFile.is_open()) {

		while(getline(listFile,modelName)) { //Retrieve the name of the model
			modelFull = modelsPath + modelName;
			//cout << "ModelFull: " << modelFull << "\n";
			modelFile.open(modelFull.c_str()); //Open the file that contains the probabilities

			if (modelFile.is_open()) {
				try { //In case there's something in the model's folder that shouldn't be there
					cout << "Model: " << modelName << "\n";
					while(modelFile >> index >> value) {
						model.push_back(value); //Store the model values
						//cout << "Model value: " << value << "\n";
					}
					//cout << "Model size: " << model.size() << "\n";

					//For each read calculate the score for the model
					for(int i=0; i<reads->getSequences().size(); i++) {
						tmpScore = this->scoreRead((string)reads->getSequences().at(i), order, kmers->getKmerList(), model);
						cout << "Score for read "<< i << ": " << tmpScore << "\n";

						//Replace the score stored if the new score is higher
						if(this->scores.size() < reads->getSequences().size()) {
							this->scores.push_back(tmpScore);
							this->modelNames.push_back((string)modelName.substr(0,modelName.find(".")));
						}
						else {
							if (tmpScore > this->scores.at(i)) {
								this->scores.at(i) = tmpScore;
								this->modelNames.at(i) = modelName.substr(0,modelName.find("."));
							}
						}

					} //End while scoring reads
				} catch(...) {}

				modelFile.close();
				model.clear();
			} //End if model was loaded
		} //End while reading models

		//Write the final scores to a file
		scoreResults.open(outputFile.c_str());
		if (scoreResults.is_open()) {
			scoreResults << "Best score\tBest model\n";
			for(int i=0; i<this->scores.size(); i++) {
				scoreResults << this->scores.at(i) << "\t" << this->modelNames.at(i) << "\n";
			}
			scoreResults.close();
		}

		listFile.close();
	} //End if list of models was read

	reads->~CReads();
	kmers->~CKmers();
}

/*
int main(int argc, char* argv[]) {

	//string pathToModels = "/scratch/giri_projects/vanessa/Azad/scripts/model_database/";
	//string pathToModels = "/Users/vanessa/Documents/Work/ResearchInProgress/BioRG/Metagenomics/smallList/";
	string pathToModels = "/Users/vanessa/Documents/Work/ResearchInProgress/BioRG/Metagenomics/signature_6order/";
	CScore* s = new CScore();
	struct timeval tv1, tv2;
	gettimeofday(&tv1, NULL);

	//for(int i=0; i<1000; i++) {
		//s->scoreModels(pathToModels,"test.fa","scores.txt",6);
	//}

	s->scoreModels(pathToModels,"test.fa","scores.txt",6);
	//s->scoreModels(pathToModels,"test.fa","scores.txt",8);

	gettimeofday(&tv2, NULL);
	double tm = (double) (tv2.tv_usec - tv1.tv_usec)/1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
	cout << "Time taken in execution = " << tm << " seconds\n";

	return 0;
}
*/
