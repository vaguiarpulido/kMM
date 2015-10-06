// TMC Faster
__constant__ int mapping[20] = {0, -1, 3, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1};

// Thread i will score genome[i*seqlength] to genome[i*seqlength+(seqlength-1)]
__global__ void scoreReads(char* genome, int seqLength, int order, float* model, float* scores) {
   int i = blockIdx.x;  // Thread identifier, assign to i
   int j = threadIdx.x;
   

   // Keep scores in shared memory
   extern __shared__ float kmer_scores[];  // Call this with [lengths[i] / order + 1];

   // Start spot
   int seqspot = i*seqLength;
   int startspot = seqspot+(j*order);

   // Quick loop, check for n's
   // Actually, decided to inline it rather than loop twice.
   int a;
   bool nFlag = false;
   int mapVal = 0;
   for (a = startspot; a < startspot+order; a++) {
      if (genome[a] == 'N') { 
         kmer_scores[j] = 0;
         nFlag = true;
         break;
      }
      else 
         mapVal = 4*mapVal + mapping[(int)genome[a]-65];
   }

   if (!nFlag)
      kmer_scores[j] = model[mapVal];

   __syncthreads();
   // Do the addition in parallel as well.
   int k = seqLength / 2;
   while (k >= 1) {
      kmer_scores[j] = kmer_scores[j] + kmer_scores[j+k];
      k /= 2;
      __syncthreads();
   }

   // The first kmer_score is now the final.
   scores[i] = kmer_scores[0];
}


