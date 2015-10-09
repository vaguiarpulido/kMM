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

#include "scoreReads.h"
#include <math.h>

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
	int index=-1;//, mappedIndex=-1;
	string tmpKmer="", tmpRead="";
	//float tmpScore=0.0;
	//Prepare to load the reads
	CSequences* reads = new CSequences(readsFileName);


        // GPU arrays
        int num_seq = reads->getSequences().size();
        int read_length = reads->getSequences()[0].size(); // TMC for now assuming same length
        int nucleotides = num_seq*read_length;
	float* cpu_scores = (float*) malloc(num_seq*sizeof(float));
 
        char* cpu_genome = (char*) malloc(nucleotides*sizeof(char));
        // data() does not work, have to copy manually for now.
        for (int i = 0; i < num_seq; i++) {
	   cout << "Sequence " << i << ": " << reads->getSequences()[i] << endl;
           for (int j = 0; j < read_length; j++) {
              cpu_genome[i*read_length+j] = reads->getSequences()[i][j];
           }
        }
        cout << "CPU GENOME: " << string(cpu_genome) << endl;

        char* gpu_genome;
        cudaMalloc((void**) &gpu_genome, nucleotides*sizeof(char));
        cudaMemcpy(gpu_genome, cpu_genome, nucleotides*sizeof(char), cudaMemcpyHostToDevice); // TMC I know this works for vectors of ints, need to check vectors of strings

	//Prepare to get the list of possible kmers for a model
	CKmers* kmers = new CKmers(order);

        order++;
	//Get the full list of models
	if(modelsPath.compare(modelsPath.length()-1,1,"/") !=0) {
		modelsPath += "/";
	}
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
					int num_models = pow(4,order) + pow(4, order-1);
                                        float* cpu_model = (float*) malloc(num_models * sizeof(float));
					cout << "Model: " << modelName << "\n";
                                        int i=0;
					while(i < num_models && modelFile >> index >> value) {
                                                cpu_model[i] = value;
						//model.push_back(value); //Store the model values
						//cout << "Model value: " << value << "\n";
                                                i++;
					}
					//cout << "Model size: " << model.size() << "\n";
					//cout << "First element: " << model.data()[0] << "\n";
					float* gpu_model;
					cudaMalloc((void**) &gpu_model, num_models*sizeof(float));
					cudaMemcpy(gpu_model, cpu_model, num_models*sizeof(float), cudaMemcpyHostToDevice);


           				float* gpu_scores;
					cudaMalloc((void**) &gpu_scores, num_seq*sizeof(float));

					//For each read calculate the score for the model
					// Call with num_seq blocks of order threads.
					int num_kmers = read_length - order + 1 + 1;
					//printf("Calling kernel.\n");
					//cout << "Using " << num_seq << " blocks of " << num_kmers << " threads.  Shared memory contains " << num_kmers << " floats." << endl;
					scoreReads<<<num_seq, num_kmers, num_kmers*sizeof(float)>>>(gpu_genome, read_length, order, gpu_model, gpu_scores);
					cudaDeviceSynchronize();
					//printf("Called kernel.\n");
					cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    printf("CUDA error: %s\n", cudaGetErrorString(error));
    exit(-1);
  }

					cudaMemcpy(cpu_scores, gpu_scores, num_seq*sizeof(float), cudaMemcpyDeviceToHost);
					//if (cpu_scores[0] != cpu_scores[1])
                                        //   cout << "Warning: Sequence 0 has " << cpu_scores[0] << " and Sequence 1 has " << cpu_scores[1] << endl;
					cudaFree(gpu_model);
					cudaFree(gpu_scores);
					free(cpu_model);
					for(int i=0; i<reads->getSequences().size(); i++) {
						//tmpScore = this->scoreRead((string)reads->getSequences().at(i), order, kmers->getKmerList(), model);
						//cout << "Score for read "<< i << ": " << cpu_scores[i] << "\n";

						//Replace the score stored if the new score is higher
						// TMC for now removing, since we are only doing one model
						//cout << "Score for sequence " << i << " (press return to continue): " << cpu_scores[i] << endl;
						//int x;cin >> x;
						if(this->scores.size() < reads->getSequences().size()) {
							this->scores.push_back(cpu_scores[i]);
							this->modelNames.push_back((string)modelName.substr(0,modelName.find(".")));
						}
						else {
							if (cpu_scores[i] > this->scores.at(i)) {
								this->scores.at(i) = cpu_scores[i];
								this->modelNames.at(i) = modelName.substr(0,modelName.find("."));
							}
						}
			//			exit(1);
					} //End while scoring reads
				} catch(...) {}

				modelFile.close();
				//cout << "Model cleared." << endl;
				model.clear();
			} //End if model was loaded
		} //End while reading models

		//Write the final scores to a file
		scoreResults.open(outputFile.c_str());
		if (scoreResults.is_open()) {
			scoreResults << "Best score\tBest model\n";
			for(int i=0; i < this->scores.size(); i++) {
				scoreResults << this->scores.at(i) << "\t" << this->modelNames.at(i) << "\n";
			}
			scoreResults.close();
		}

		listFile.close();
	} //End if list of models was read

	reads->~CSequences();
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
