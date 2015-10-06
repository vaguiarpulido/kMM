// TMC Faster
__constant__ int mapping[20] = {0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

__device__ int mapKmer(char * str, int loc, int len) {
  // Takes a substring of str starting from loc of length len and computes its mapped value
  int mapVal = 0;
  for(int i = loc; i < len; i++) {
    mapVal = 4 * mapVal + mapping[(int)str[i]-65];
  }
  return mapVal;
}

/*__global__ void removeNs(char* genome, int seqLength, char** sequences, int* lengths) {
   int i = threadIdx.x;

   int startspot = i*seqLength;
   int endspot = startspot + seqLength - 1;
 
   int j;
   int count = 0;
   for (j = startspot; j <= endspot; j++) {
      if (genome[j] != 'N') {
         sequences[i][count] = genome[j];
         count++;
      }
   }
   lengths[i] = count;
}*/

// Thread i will score genome[i*seqlength] to genome[i*seqlength+(seqlength-1)]
__global__ void scoreReads(char** sequences, int* lengths, int order, float* model, float* scores) {
   int i = blockIdx.x;  // Thread identifier, assign to i
   int j = threadIdx.x;
   

   // Keep scores in shared memory
   extern __shared__ float kmer_scores[];  // Call this with [lengths[i] / order + 1];

   // Start spot
   int startspot = j*order;

   // Quick loop, check for n's
   int a;
   bool flag = false;
   for (a = startspot; a < startspot+order; a++) {
      if (sequences[i][a] == 'N') { 
         flag = true;
         break;
      }
   }

   if (flag)
      kmer_scores[j] = 0;
   else
      kmer_scores[j] = model[mapKmer(sequences[i], j*order, order)];

   __syncthreads();
   // Do the addition in parallel as well.
   int k = lengths[i] / 2;
   while (k >= 1) {
      kmer_scores[j] = kmer_scores[j] + kmer_scores[j+k];
      k /= 2;
      __syncthreads();
   }

   // The first kmer_score is now the final.
   scores[i] = kmer_scores[0];
}


