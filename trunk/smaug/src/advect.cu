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
#include "../include/gradops_advect.cuh"
#include "../include/dervfields_advect.cuh"



__global__ void advect_parallel(struct params *p,  real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, int order)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;
  int order1;
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
  int ix[NDIM];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  

   j=iindex/ni;
   //i=iindex-j*(iindex/ni);
   i=iindex-(j*ni);
   ix[0]=i;
   ix[1]=j;
   if(order==0 || order==1)
     dt=(p->dt)/2.0;
   if(order==2 )
     dt=(p->dt);
   if(order==3 )
     dt=(p->dt)/6.0;

   if(order==3)
     order1=0;
   else
     order1=order+1;

  //order1=1;
//dt=(p->dt);
  if(p->rkon != 1)
  {
    dt=(p->dt);
    order1=0;
  }


  //advance the solution for one of the advect steps
  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{ 
   
		for(int f=rho; f<=b3; f++)           
 			//wmod[fencode_advect(p,i,j,f)]=((w[fencode_advect(p,i+1,j,f)]+w[fencode_advect(p,i-1,j,f)]+w[fencode_advect(p,i,j+1,f)]+w[fencode_advect(p,i,j-1,f)])/4.0)+dt*dwn1[(NVAR*ni*nj*(order-1))+fencode_advect(p,i,j,f)];
                   wmod[NVAR*(ni*nj*(order1))+fencode_advect(p,i,j,f)]=wmod[NVAR*(ni*nj*(order))+fencode_advect(p,i,j,f)]+dt*dwn1[fencode_advect(p,i,j,f)];
	}

 __syncthreads();



  if((p->rkon == 1) && i<((p->n[0])) && j<((p->n[1])))
	{ 
   
		for(int f=rho; f<=b3; f++) 
                {   
                       if(order==1)
			   wmod[fencode_advect(p,i,j,f)]=wmod[NVAR*(ni*nj*(order1))+fencode_advect(p,i,j,f)];	
                        /*else if(order==1) 
			   wmod[NVAR*(ni*nj*(2))+fencode_advect(p,i,j,f)]=w[fencode_advect(p,i,j,f)];*/	
                        if(order==2)  
			   wmod[NVAR*(ni*nj*(1))+fencode_advect(p,i,j,f)]=(wmod[NVAR*(ni*nj*(1))+fencode_advect(p,i,j,f)]	+2*wmod[NVAR*(ni*nj*(2))+fencode_advect(p,i,j,f)]    + wmod[NVAR*(ni*nj*(3))+fencode_advect(p,i,j,f)]- 4*w[fencode_advect(p,i,j,f)])/3;

                        else if(order==3)      
 			//wmod[fencode_advect(p,i,j,f)]=((w[fencode_advect(p,i+1,j,f)]+w[fencode_advect(p,i-1,j,f)]+w[fencode_advect(p,i,j+1,f)]+w[fencode_advect(p,i,j-1,f)])/4.0)+dt*dwn1[(NVAR*ni*nj*(order-1))+fencode_advect(p,i,j,f)];
                   wmod[fencode_advect(p,i,j,f)]=wmod[NVAR*(ni*nj)+fencode_advect(p,i,j,f)]+wmod[fencode_advect(p,i,j,f)];
                }
	}

 __syncthreads();



  if((order==3) && (p->rkon == 1) && i<((p->n[0])) && j<((p->n[1])))   
		for(int f=rho; f<=b3; f++) 
                   wmod[fencode_advect(p,i,j,f)]=wmod[NVAR*(ni*nj)+fencode_advect(p,i,j,f)]+wmod[fencode_advect(p,i,j,f)];



 __syncthreads();

if(i<((p->n[0])) && j<((p->n[1])))
	{		
               //for(int f=rho; f<=b3; f++)
               //{               
               //   wmod[fencode_advect(p,i,j,f)]=w[fencode_advect(p,i,j,f)];
               //   wnew[fencode_advect(p,i,j,f)]=0.0;
               //}
               for(int f=current1; f<=f3; f++)
                  wd[fencode_advect(p,i,j,f)]=0; 
        }
               __syncthreads();


 /* if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{		               
               computej_advect(wmod+((p->n[0])*(p->n[1])*order1),wd,p,i,j);
               computepk_advect(wmod+((p->n[0])*(p->n[1])*order1),wd,p,i,j);
               computept_advect(wmod+((p->n[0])*(p->n[1])*order1),wd,p,i,j);

               computebdotv_advect(wmod+((p->n[0])*(p->n[1])*order1),wd,p,i,j);
               computedivb_advect(wmod+((p->n[0])*(p->n[1])*order1),wd,p,i,j);
         }
              __syncthreads();
  if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               computec_advect(wmod+((p->n[0])*(p->n[1])*order1),wd,p,i,j);
        }
              __syncthreads();*/


  
}

/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_advect(char *label)
{
  // we need to synchronise first to catch errors due to
  // asynchroneous operations that would otherwise
  // potentially go unnoticed

  cudaError_t err;

  err = cudaThreadSynchronize();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }
}




int cuadvect(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, int order)
{


//printf("calling propagate solution\n");

    //dim3 dimBlock(blocksize, blocksize);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimproduct_advect(*p)+numThreadsPerBlock-1) / numThreadsPerBlock;

//__global__ void prop_parallel(struct params *p, real *b, real *w, real *wnew, real *wmod, 
  //  real *dwn1, real *dwn2, real *dwn3, real *dwn4, real *wd)
     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
     advect_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wnew, *d_wmod, *d_dwn1,  *d_wd, order);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called boundary\n");  
     //cudaThreadSynchronize();
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_b,*d_w,*d_wnew);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
 

  //  cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

     //following used for testing to check current soundspeeds etc
     //cudaMemcpy(*w, *d_wd, 7*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}






