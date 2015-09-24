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

CScore::CScore() {


}

CScore::~CScore() {

}

float CScore::scoreRead(string read, int order, map<string, int> kmerMap, vector<float> model) {
	float score = 0.0;
	string tmpKmer = read.substr(0, order+1); //From the beginning, take k+1 characters
	int mappedIndex = kmerMap.at(tmpKmer); //This will return the position of the kmer
	score = score + model.at(mappedIndex);

	for(int j=order; j<read.length()-order+1; j++) { //Calculate the score of a read
		//cout << "tmpKmer " << j << ": " << tmpKmer << "\n";
		tmpKmer.erase(tmpKmer.begin());
		tmpKmer.push_back((char)read.at(j));
		mappedIndex = kmerMap.at(tmpKmer); //This will return the position of the kmer
		score = score + model.at(mappedIndex);
		//cout << "Partial score "<< j << ": " << tmpScore << "\n";
	}
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

	if (listFile.is_open()) {

		while(getline(listFile,modelName)) { //Retrieve the name of the model
			modelFull = modelsPath + modelName;
			modelFile.open(modelFull.c_str()); //Open the file that contains the probabilities

			if (modelFile.is_open()) {
				//cout << "Model: " << modelName << "\n";
				while(modelFile >> index >> value) {
					model.push_back(value); //Store the model values
					//cout << "Model value: " << value << "\n";
				}

				//For each read calculate the score for the model
				for(int i=0; i<reads->getSequences().size(); i++) {
					tmpScore = this->scoreRead((string)reads->getSequences().at(i), order, kmers->getKmerList(), model);

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

				} //End for scoring reads

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


int main(int argc, char* argv[]) {

	string pathToModels = "/Users/vanessa/Documents/Work/ResearchInProgress/BioRG/Metagenomics/smallList/";
	CScore* s = new CScore();
	s->scoreModels(pathToModels,"test.txt","scores.txt",4);

	return 0;
}

