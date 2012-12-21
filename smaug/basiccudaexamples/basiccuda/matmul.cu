#include <stdio.h>

// Thread block size
#define BLOCK_SIZE 16
// Forward declaration of the device multiplication function
__global__ void Muld(float*, float*, int, int, float*, float *);




// Host multiplication function
// Compute C = A * B
//   hA is the height of A
//   wA is the width of A
//   wB is the width of B
void Mul(const float* A, const float* B, int hA, int wA, int wB,
         float* C, float *par)
{
    


    int size;
    // Load A and B to the device
    float* Ad;
    float* pard;
    size = hA * wA * sizeof(float);
    cudaMalloc((void**)&Ad, size);
    cudaMalloc((void **)&pard,2);
    cudaMemcpy(Ad, A, size, cudaMemcpyHostToDevice);
    float* Bd;
    size = wA * wB * sizeof(float);
    cudaMalloc((void**)&Bd, size);
    cudaMemcpy(Bd, B, size, cudaMemcpyHostToDevice);
    // Allocate C on the device
    float* Cd;
    size = hA * wB * sizeof(float);
    cudaMalloc((void**)&Cd, size);
    // Compute the execution configuration assuming
    // the matrix dimensions are multiples of BLOCK_SIZE
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid(wB / dimBlock.x, hA / dimBlock.y);
    // Launch the device computation
    Muld<<<dimGrid, dimBlock>>>(Ad, Bd, wA, wB, Cd,pard);
    // Read C from the device
    cudaMemcpy(C, Cd, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(par,pard,2,cudaMemcpyDeviceToHost);
    // Free device memory
    cudaFree(Ad);
    cudaFree(Bd);
    cudaFree(Cd);
}

// Device multiplication function called by Mul()
// Compute C = A * B
//   wA is the width of A
//   wB is the width of B
__global__ void Muld(float* A, float* B, int wA, int wB, float* C, float* pard)
{
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    // Index of the first sub-matrix of A processed by the block
    int aBegin = wA * BLOCK_SIZE * by;
    // Index of the last sub-matrix of A processed by the block
    int aEnd   = aBegin + wA - 1;
    // Step size used to iterate through the sub-matrices of A
    int aStep = BLOCK_SIZE;
    // Index of the first sub-matrix of B processed by the block
    int bBegin = BLOCK_SIZE * bx;
    // Step size used to iterate through the sub-matrices of B
    int bStep = BLOCK_SIZE * wB;
    // The element of the block sub-matrix that is computed
    // by the thread
    float Csub = 0;
     if(bx<1 ) pard[0]=bx;
     if(tx<1) pard[1]=tx;
    // Loop over all the sub-matrices of A and B required to
    // compute the block sub-matrix
    for (int a = aBegin, b = bBegin;
             a <= aEnd;
             a += aStep, b += bStep) {
        // Shared memory for the sub-matrix of A
        __shared__ float As[BLOCK_SIZE][BLOCK_SIZE];
        // Shared memory for the sub-matrix of B
        __shared__ float Bs[BLOCK_SIZE][BLOCK_SIZE];
        // Load the matrices from global memory to shared memory;
        // each thread loads one element of each matrix
        As[ty][tx] = A[a + wA * ty + tx];
        Bs[ty][tx] = B[b + wB * ty + tx];
        // Synchronize to make sure the matrices are loaded
        __syncthreads();
        // Multiply the two matrices together;
        // each thread computes one element
        // of the block sub-matrix
        for (int k = 0; k < BLOCK_SIZE; ++k)
          Csub += As[ty][k] * Bs[k][tx];
      // Synchronize to make sure that the preceding
      // computation is done before loading two new
      // sub-matrices of A and B in the next iteration
      __syncthreads();
  }
  // Write the block sub-matrix to global memory;
  // each thread writes one element
  int c = wB * BLOCK_SIZE * by + BLOCK_SIZE * bx;
  C[c + wB * ty + tx] = Csub;
}


int main(void)
{
    int hA,wA,wB;
    int size;
    int i;

    hA=20;
    wA=4;
    wB=20;
    // Load A and B to the device
    float* A;
    float *par;
    size = hA * wA * sizeof(float);
    A=(float *)malloc( size);
    for(i=0; i<hA*wA; i++) A[i]=i;
    
    float* B;
    size = wA * wB * sizeof(float);
    B=(float *)malloc( size);
    for(i=0; i<wA*wB; i++) B[i]=2*i;
    // Allocate C on the device
    float* C;
    size = hA * wB * sizeof(float);
    C=(float *)malloc(size);


    par=(float *)malloc(1);




    Mul(A,B,hA,wA,wB,C,par);
    for(i=0; i<hA*wA; i++) printf("%d %f ",i,A[i]);
    printf("\n");
    for(i=0; i<wA*wB; i++) printf("%d %f ",i,B[i]);
    printf("\n");
    for(i=0; i<hA*wB; i++) printf("%d %f ",i,C[i]);
    printf("\n");



    printf("\n%f %f\n",par[0],par[1]);

}

