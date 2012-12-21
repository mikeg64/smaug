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
 #ifdef USE_SAC_3D
   int kp,kpg;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

 


int shift=order*NVAR*dimp;

  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif    


     #ifdef USE_SAC_3D
           for( f=rho; f<=b3; f++)
     #else
           for( f=rho; f<=b2; f++)
     #endif
             {  
         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{

 

              // bc3_periodic1_b(wmod+order*NVAR*dimp,p,iia,f);  //for OZT
              
     #ifdef USE_SAC_3D
         ;//if((f!=mom1 || f !=mom2 || f != mom3) && (p->it)>0)      
      #else
         ;//if(f!=mom1 || f !=mom2 )
      #endif             
                // bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
                

                //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);

                  bc3_fixed_dir_b(wmod+order*NVAR*dimp,p,bp,iia,f,dir);

	}

               }

 __syncthreads();
 

}



__global__ void setboundary_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  
                real val=0;
  int nk;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
   nk=p->n[2];
#endif  
  

   /* real *wtest;
   int i1,i2,i3;

    switch(dir)
    {
     case 0:
	i1=127;
	i2=63;
	i3=63;
     break;
     case 1:
	i1=63;
	i2=127;
	i3=63;
     break;

     case 2:
	i1=63;
	i2=63;
	i3=127;
     break;


    }*/

  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif     

int shift=order*NVAR*dimp;


         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{

 

              // bc3_periodic1_b(wmod+order*NVAR*dimp,p,iia,f);  //for OZT
              
     #ifdef USE_SAC_3D
         ;//if((f!=mom1 || f !=mom2 || f != mom3) && (p->it)>0)      
      #else
         ;//if(f!=mom1 || f !=mom2 )
      #endif             

		  bc3_setfixed_dir_b(wmod,p,bp,iia,field,dir);

                  /*if(i==i1 && j==i2 && k==i3)
                  {
                    wtest=wmod;//+order*NVAR*dimp;
                    p->test=wmod[encode3_b(p,i,j,k,field)];

                   }*/
	}

 __syncthreads();
 

}





int cuboundary(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp, struct state **d_s,  real **d_wmod,  int order,int idir,int field)
{


 dim3 dimBlock(dimblock, 1);

int numBlocks = ((dimproduct_b(*p)+numThreadsPerBlock-1)) / numThreadsPerBlock;
  



int i1,i2,i3;



if(((*p)->it)==-1)
{
    //cudaMemcpy(*d_p,*p, sizeof(struct params), cudaMemcpyHostToDevice);
    setboundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);

    //cudaMemcpy(*bp,*d_bp, sizeof(struct bparams), cudaMemcpyDeviceToHost);
    //cudaMemcpy(*p,*d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

     //
    /*switch(idir)
    {
     case 0:
	i1=127;
	i2=63;
	i3=63;

      printf("dir=%d field=%d value=%20.10g  fixed=%20.10g\n",idir,field,(*p)->test,(*bp)->fixed1[encodefixed13_b(*p,1+(*p)->n[0]-i1,i2,i3,field)]);
     break;
     case 1:

	i1=63;
	i2=127;
	i3=63;

      printf("dir=%d field=%d value=%20.10g  fixed=%20.10g\n",idir,field,(*p)->test,(*bp)->fixed2[encodefixed23_b(*p,i1,1+(*p)->n[1]-i2,i3,field)]);
     break;
     case 2:

	i1=63;
	i2=63;
	i3=127;

      printf("dir=%d field=%d value=%20.10g  fixed=%20.10g\n",idir,field,(*p)->test,(*bp)->fixed3[encodefixed33_b(*p,i1,i2,1+(*p)->n[2]-i3,field)]);
     break;

    }*/


}
else
{
    boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);

    cudaThreadSynchronize();
    boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
    cudaThreadSynchronize();
#ifdef USE_SAC_3D
    boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
    cudaThreadSynchronize();
#endif
}
    


}








