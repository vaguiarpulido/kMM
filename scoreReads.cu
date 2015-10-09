#include <stdio.h>

// TMC Faster
__constant__ int mapping[20] = {0, -1, 3, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1};

// Thread i will score genome[i*seqlength] to genome[i*seqlength+(seqlength-1)]
__global__ void scoreReads(char* genome, int seqLength, int order, float* model, float* scores) {
   int i = blockIdx.x;  // Thread identifier, assign to i
   int j = threadIdx.x;
   // Keep scores in shared memory
   extern __shared__ float kmer_scores[];  // Call this with [lengths[i] / order + 1];

   //if (i ==0) printf("%s\n", genome);
   // Start spot
   int seqspot = i*seqLength;
   int startspot, stopspot;
   if (j == 0) {
      startspot = seqspot;
      stopspot = startspot+(order-1);
   }
   else{
      startspot = seqspot+(j-1);
      stopspot = startspot+order;
   }
   //printf("Block %d Thread %d Startspot %d Stopspot %d\n", i, j, startspot, stopspot);
   // Quick loop, check for n's
   // Actually, decided to inline it rather than loop twice.
   int a;
   bool nFlag = false;
   int mapVal = 0;
   for (a = startspot; a < stopspot; a++) {
      //if (j == 0 && i == 0) printf("%d %d\n", startspot, a);
      //if (j == 0 && i == 0) printf("%c\n", genome[a]);
      if (genome[a] == 'N') { 
         //if (i == 0) printf("FOUND N, BREAKING.\n");
         kmer_scores[j] = 0;
         nFlag = true;
         break;
      }
      else 
         mapVal = 4*mapVal + mapping[(int)genome[a]-65];
   }
   if (j == 1 && i == 0) printf("\n");
   if (!nFlag) {
      if (j == 0) {mapVal += pow(4.0, 9.0);}
      kmer_scores[j] = model[mapVal]; // Illegal here
      //if (i == 0) printf("Thread: %d  Mapval: %d  Score: %f\n", j, mapVal, kmer_scores[j]);
    }
   __syncthreads();

    /////// TMC TAKE OUT LATER
   /*if (j == 0) {
      int m;
      float tmpscore=0;
      for (m = seqspot; m < seqspot+seqLength-order+1+1; m++) {
         if (i == 0) printf("%d: Score: %f  New Partial Score: %f\n", m+8, kmer_scores[m-seqspot], tmpscore); 
         tmpscore += kmer_scores[m-seqspot];
      }
      if (i == 0) printf("The score for sequence %d should be: %f\n", i, tmpscore);
     
   }
   __syncthreads();*/
   //////////////////////////


   // Do the addition in parallel as well.
   //if (j == 0 && i == 0) printf("Number: %d", seqLength-order+1);
   if ((j == 0) && (((seqLength-order+1+1) % 2 == 1))) {/*if (i == 0) printf("Thread %d adding %d (%f) and %d (%f)", j, j, kmer_scores[j], seqLength-order+1, kmer_scores[seqLength-order+1]);*/ kmer_scores[j] = kmer_scores[j] + kmer_scores[seqLength-order+1]; /*if (i == 0) printf(" to get: %f\n", kmer_scores[j]);*/ }
   int k = (seqLength-order+1+1)/2;
   while (k >= 1) {
      //if (i == 0 && j == 0) printf("k is %d\n", k);
      if (j < k) {
	 //if (i == 0) printf("k: %d Thread %d adding %d (%f) and %d (%f) to get: %f\n", k, j, j, kmer_scores[j], j+k, kmer_scores[j+k], kmer_scores[j]+kmer_scores[j+k]);
         kmer_scores[j] = kmer_scores[j] + kmer_scores[j+k]; // Illegal here
      }
      __syncthreads();
      if (k != 1 && k % 2 != 0 && j == 0) {/*printf("k: %d Thread %d adding %d (%f) and %d (%f) to get: %f\n", k, j, j, kmer_scores[j], k-1, kmer_scores[k-1], kmer_scores[j]+kmer_scores[k-1]);*/              kmer_scores[j] = kmer_scores[j] + kmer_scores[k-1]; }// Uneven split, would be left out otherwise.
      k /= 2;
   }

   // The first kmer_score is now the final.
   if (j == 0) {
      scores[i] = kmer_scores[0];
      //printf("Kernel score for sequence %d: %f\n", i, scores[i]);
   }
}


