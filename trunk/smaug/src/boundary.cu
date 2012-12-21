#include "../include/cudapars.h"
#include "../include/paramssteeringtest1.h"

/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "../include/smaugcukernels.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "../include/gradops_b.cuh"

__global__ void boundary_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  
  int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));


 #if defined USE_SAC
    j=iindex/(ni);
   i=iindex-(jp*ni);
#endif  


int shift=order*NVAR*dimp;



     iia[0]=i;
     iia[1]=j;
     k=0;

      for( f=rho; f<=b2; f++)
      {  
	       if(i<((p->n[0])) && j<((p->n[1]))) 
		{
		       bc3_periodic1_dir_b(wmod+order*NVAR*dimp,p,iia,f,dir);  //for OZT
		       //  bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
		       //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);
		}
       }

 __syncthreads();





  
}

int cuboundary(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp, struct state **d_s,  real **d_wmod,  int order,int idir,int field)
{


 dim3 dimBlock(dimblock, 1);

int numBlocks = ((dimproduct_b(*p)+numThreadsPerBlock-1)) / numThreadsPerBlock;
if((*p)->boundtype[field][0]==0)
    boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);

    cudaThreadSynchronize();
if((*p)->boundtype[field][1]==0)
    boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
    cudaThreadSynchronize();
}
