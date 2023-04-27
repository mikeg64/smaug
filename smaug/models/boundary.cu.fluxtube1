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

__global__ void boundary0_contcd4_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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

f=field;
  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif


     #ifdef USE_SAC_3D
          // for( f=rho; f<=b3; f++)
          //   {
          k=iia[2];
                if((i==0 || i==1))
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,4-i,j,k,f)];
                else if((( i==((p->n[0])-1)   )))
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i-4,j,k,f)];
                else if(((  i==((p->n[0])-2) )) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i-2,j,k,f)];
          //          }

    #else
          // for( f=rho; f<=b2; f++)
          //   {

                if((i==0 || i==1) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,4-i,j,k,f)];
                else if((( i==((p->n[0])-1)   )) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i-4,j,k,f)];
                else if(((  i==((p->n[0])-2) ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i-2,j,k,f)];
		//}


    #endif

 __syncthreads();


}



__global__ void boundary0_cont_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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

f=field;
  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif







     #ifdef USE_SAC_3D
          // for( f=rho; f<=b3; f++)
          //   {
          k=iia[2];
                if((i==0 || i==1))
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,3,j,k,f)];
                else if((( i==((p->n[0])-1)   )))
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,(p->n[0])-3,j,k,f)];
                else if(((  i==((p->n[0])-2) )) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,(p->n[0])-3,j,k,f)];
          //          }

    #else
          // for( f=rho; f<=b2; f++)
          //   {

                if((i==0 || i==1) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,3,j,k,f)];
                else if((( i==((p->n[0])-1)   )) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,(p->n[0])-3,j,k,f)];
                else if(((  i==((p->n[0])-2) ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,(p->n[0])-3,j,k,f)];
		//}


    #endif

 __syncthreads();


}




__global__ void boundary1_contcd4_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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
f=field;

  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif


     #ifdef USE_SAC_3D
        //   for( f=rho; f<=b3; f++)
        //     {
          k=iia[2];
                if((j==0 || j==1) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,4-j,k,f)];
                else if((( j==((p->n[1])-1)   ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j-4,k,f)];
                else if(((  j==((p->n[1])-2) ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j-2,k,f)];

            //        }

    #else
         //  for( f=rho; f<=b2; f++)
         //    {

                if((j==0 || j==1) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,4-j,k,f)];
                else if((( j==((p->n[1])-1)   ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j-4,k,f)];
                else if(((  j==((p->n[1])-2) ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j-2,k,f)];
		//}


    #endif

 __syncthreads();


}





__global__ void boundary1_cont_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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
f=field;

  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif


     #ifdef USE_SAC_3D
        //   for( f=rho; f<=b3; f++)
        //     {
          k=iia[2];
                if((j==0 || j==1) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,3,k,f)];
                else if((( j==((p->n[1])-1)   ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,((p->n[1])-3),k,f)];
                else if(((  j==((p->n[1])-2) ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,((p->n[1])-3),k,f)];

            //        }

    #else
         //  for( f=rho; f<=b2; f++)
         //    {

                if((j==0 || j==1) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,3,k,f)];
                else if((( j==((p->n[1])-1)   ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,(p->n[1])-2,k,f)];
                else if(((  j==((p->n[1])-2) ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,(p->n[1])-2,k,f)];
		//}


    #endif

 __syncthreads();


}








__global__ void boundary2_contcd4_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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
  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif
f=field;


     #ifdef USE_SAC_3D
         //  for( f=rho; f<=b3; f++)
         //    {
          k=iia[2];
               if((k==0 || k==1) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j,4-k,f)];
                else if((( k==((p->n[2])-1)   ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j,k-4,f)];
                else if(((  k==((p->n[2])-2) ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j,k-2,f)];

           //         }



    #endif

 __syncthreads();


}












__global__ void boundary2_cont_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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
  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif
f=field;


     #ifdef USE_SAC_3D
         //  for( f=rho; f<=b3; f++)
         //    {
          k=iia[2];
               if((k==0 || k==1) )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j,3,f)];
                else if((( k==((p->n[2])-1)   ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j,(p->n[2])-3,f)];
                else if(((  k==((p->n[2])-2) ))  )
                    wmod[  shift+encode3_b(p,i,j,k,f)]=wmod[  shift+encode3_b(p,i,j,(p->n[2])-3,f)];

           //         }



    #endif

 __syncthreads();


}








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
f=field;

  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif


     /*#ifdef USE_SAC_3D
           for( f=rho; f<=b3; f++)
     #else
           for( f=rho; f<=b2; f++)
     #endif
             { */
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

          //     }

 __syncthreads();


}





__global__ void boundary0_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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
f=field;

  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif


    // #ifdef USE_SAC_3D
      //     for( f=rho; f<=b3; f++)
     //#else
     //      for( f=rho; f<=b2; f++)
     //#endif
       //      {
		        if(i==0 || i==1 )
		          wmod[encode3_b(p,i,j,k,f)+order*NVAR*dimp]=bp->fixed1[encodefixed13_b(p,i,j,k,f)];
		        if((i==((p->n[0])-1)) || (i==((p->n[0])-2)) )
		          wmod[encode3_b(p,i,j,k,f)+order*NVAR*dimp]=bp->fixed1[encodefixed13_b(p,1+(p->n[0])-i,j,k,f)];
	 //}

 __syncthreads();


}


__global__ void boundary1_per_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];

  int ip,jp,kp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);

   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif


f=field;

int shift=order*NVAR*dimp;



   //  iia[0]=ip;
  //   iia[1]=jp;
   //  k=0;

      //for( f=rho; f<=b2; f++)
      //{

	      // if(i<((p->n[0])) && j<((p->n[1])))
		//{
		       //bc3_periodic1_dir2_b(wmod+order*NVAR*dimp,p,iia,f,dir);  //for OZT
		       //  bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
		       //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);
		//}
       //}



f=field;
//kp=0;


                if((jp==0 || jp==1)  && ip>1 && ip<((p->n[1])-2) )
                  wmod[shift+encode3_b(p,ip,jp,kp,f)]=wmod[shift+encode3_b(p,ip,(p->n[1])-4+jp,kp,f)];

                else if(((jp==((p->n[1])-1)) || (jp==((p->n[1])-2))) && ip>1 && ip<((p->n[1])-2) )
                  wmod[shift+encode3_b(p,ip,jp,kp,f)]=wmod[shift+encode3_b(p,ip,4-(p->n[1])+jp,kp,f)];




 __syncthreads();






}

__global__ void boundary2_per_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];

  int ip,jp,kp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);

   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif



int shift=order*NVAR*dimp;



   //  iia[0]=ip;
  //   iia[1]=jp;
   //  k=0;

      //for( f=rho; f<=b2; f++)
      //{

	      // if(i<((p->n[0])) && j<((p->n[1])))
		//{
		       //bc3_periodic1_dir2_b(wmod+order*NVAR*dimp,p,iia,f,dir);  //for OZT
		       //  bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
		       //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);
		//}
       //}



f=field;
//kp=0;


                if((kp==0 || kp==1) && ip>1 && ip<((p->n[1])-2) )
                  wmod[shift+encode3_b(p,ip,jp,kp,f)]=wmod[shift+encode3_b(p,ip,jp,(p->n[2])-4+kp,f)];

                else if(((kp==((p->n[2])-1)) || (kp==((p->n[2])-2))) && ip>1 && ip<((p->n[1])-2))
                  wmod[shift+encode3_b(p,ip,jp,kp,f)]=wmod[shift+encode3_b(p,ip,jp,4-(p->n[2])+kp,f)];




 __syncthreads();






}

__global__ void boundary1_fix_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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
  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif
f=field;


    // #ifdef USE_SAC_3D
      //     for( f=rho; f<=b3; f++)
     //#else
     //      for( f=rho; f<=b2; f++)
     //#endif
     //        {
		        if(j==0 || j==1 )
		          wmod[encode3_b(p,i,j,k,f)+order*NVAR*dimp]=bp->fixed1[encodefixed23_b(p,i,j,k,f)];
		        if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )
		          wmod[encode3_b(p,i,j,k,f)+order*NVAR*dimp]=bp->fixed1[encodefixed23_b(p,i,1+(p->n[1])-j,k,f)];
	// }

 __syncthreads();


}



__global__ void boundary1_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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


f=field;


int shift=order*NVAR*dimp;
  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif


   //  #ifdef USE_SAC_3D
   //        for( f=rho; f<=b3; f++)
   //          {
		        if(j==0 || j==1 )
		          wmod[encode3_b(p,i,j,k,f)+order*NVAR*dimp]=bp->fixed1[encodefixed13_b(p,i,j,k,f)];
		        if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )
		          wmod[encode3_b(p,i,j,k,f)+order*NVAR*dimp]=bp->fixed1[encodefixed13_b(p,i,1+(p->n[1])-j,k,f)];
	// }
    //#endif

 __syncthreads();


}




__global__ void boundary2_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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


f=field;


int shift=order*NVAR*dimp;
  k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif


    // #ifdef USE_SAC_3D
    //       for( f=rho; f<=b3; f++)
    //         {
		        if(k==0 || k==1 )
		          wmod[encode3_b(p,i,j,k,f)+order*NVAR*dimp]=bp->fixed1[encodefixed33_b(p,i,j,k,f)];
		        if((k==((p->n[2])-1)) || (k==((p->n[2])-2)) )
		          wmod[encode3_b(p,i,j,k,f)+order*NVAR*dimp]=bp->fixed1[encodefixed33_b(p,i,j,1+(p->n[2])-k,f)];
	// }
    //#endif

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


  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif
f=field;

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


__global__ void setboundary0_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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


k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif

int shift=order*NVAR*dimp;
f=field;


  //   #ifdef USE_SAC_3D
  //         for( f=rho; f<=b3; f++)
  //   #else
  //         for( f=rho; f<=b2; f++)
  //   #endif
  //           {
                if(i==0 || i==1 )
                  bp->fixed1[encodefixed13_b(p,i,j,k,f)]=wmod[encode3_b(p,i,j,k,f)];
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2))  )
                  bp->fixed1[encodefixed13_b(p,1+(p->n[0])-i,j,k,f)]=wmod[encode3_b(p,i,j,k,f)];

//		}
		 __syncthreads();


}



__global__ void setboundary1_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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
  f=field;


k=0;
  #ifdef USE_SAC_3D
   iia[2]=k=iindex/(nj*ni);
   iia[1]=j=(iindex-(k*(nj*ni)))/ni;
   iia[0]=i=iindex-(k*nj*ni)-(j*ni);

#else
    iia[1]=j=iindex/ni;
   iia[0]=i=iindex-(j*ni);
#endif

int shift=order*NVAR*dimp;


     //#ifdef USE_SAC_3D
     //      for( f=rho; f<=b3; f++)
     //#else
     //      for( f=rho; f<=b2; f++)
     //#endif
     //        {
                if(j==0 || j==1 )
                  bp->fixed1[encodefixed23_b(p,i,j,k,f)]=wmod[encode3_b(p,i,j,k,f)];
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2))  )
                  bp->fixed1[encodefixed23_b(p,i,1+(p->n[1])-j,k,f)]=wmod[encode3_b(p,i,j,k,f)];

	//	}
		 __syncthreads();


}


__global__ void setboundary2_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
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

f=field;

k=0;
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
       //    for( f=rho; f<=b3; f++)
       //      {
                if(k==0 || k==1 )
                  bp->fixed1[encodefixed33_b(p,i,j,k,f)]=wmod[encode3_b(p,i,j,k,f)];
                else if((k==((p->n[2])-1)) || (k==((p->n[2])-2))  )
                  bp->fixed1[encodefixed33_b(p,i,j,1+(p->n[2])-k,f)]=wmod[encode3_b(p,i,j,k,f)];

	//	}
     #endif
		 __syncthreads();


}


int cuboundary(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp, struct state **d_s,  real **d_wmod,  int order,int idir,int field)
{


 dim3 dimBlock(dimblock, 1);

int numBlocks = ((dimproduct_b(*p)+numThreadsPerBlock-1)) / numThreadsPerBlock;




int i1,i2,i3;



if(((*p)->it)==-1    )
{
    //cudaMemcpy(*d_p,*p, sizeof(struct params), cudaMemcpyHostToDevice);



	if((*p)->boundtype[field][idir][0]==4)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==0)
	boundary0_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
	    cudaThreadSynchronize();

if(idir==1)
	boundary1_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	 //   boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
            boundary2_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	   // boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	  //  cudaThreadSynchronize();
	#endif
	}






	if((*p)->boundtype[field][idir][0]==3)
	{


if(idir==0)
    setboundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);
    //cudaThreadSynchronize();
if(idir==1)
    setboundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);
    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
	    setboundary2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);
	    cudaThreadSynchronize();
	#endif
}

}
else
{

/*
	    boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
            boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);

	    cudaThreadSynchronize();
	#endif    */


	if((*p)->boundtype[field][idir][0]==5)
	{
if(idir==0)
	    boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
  //cudaThreadSynchronize();
if(idir==1)
           boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	//boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    cudaThreadSynchronize();
	    //boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	//    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
            //boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
if(idir==2)
	    boundary2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    cudaThreadSynchronize();
	#endif
	}

	if((*p)->boundtype[field][idir][0]==0)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==1)
	boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    cudaThreadSynchronize();
	    //boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    //cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
            boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    //boundary2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    cudaThreadSynchronize();
	#endif
	}


	if((*p)->boundtype[field][idir][0]==4)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==0)
	boundary0_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
	    cudaThreadSynchronize();

if(idir==1)
	boundary1_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	 //   boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
            boundary2_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	   // boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	  //  cudaThreadSynchronize();
	#endif
	}

}



}



int cuboundary1(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp, struct state **d_s,  real **d_wmod,  int order,int idir,int field)
{


 dim3 dimBlock(dimblock, 1);

int numBlocks = ((dimproduct_b(*p)+numThreadsPerBlock-1)) / numThreadsPerBlock;




int i1,i2,i3;



if(((*p)->it)==-1    )
{
    //cudaMemcpy(*d_p,*p, sizeof(struct params), cudaMemcpyHostToDevice);



	if((*p)->boundtype[field][idir][0]==4)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==0)
	;//boundary0_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
	;//    cudaThreadSynchronize();

if(idir==1)
	;//boundary1_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	 //   boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
          ;//  boundary2_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	   // boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	  //  cudaThreadSynchronize();
	#endif
	}


	if((*p)->boundtype[field][idir][0]==0)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==0)
	;//boundary0_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
	;//    cudaThreadSynchronize();

if(idir==1)
	//boundary1_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	    boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
          //  boundary2_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    cudaThreadSynchronize();
	#endif
	}




	if((*p)->boundtype[field][idir][0]==3)
	{


if(idir==0)
    setboundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);
    //cudaThreadSynchronize();
if(idir==1)
    setboundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);
    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
	    setboundary2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);
	    cudaThreadSynchronize();
	#endif
}

}
else
{

/*
	    boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
            boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);

	    cudaThreadSynchronize();
	#endif    */


	if((*p)->boundtype[field][idir][0]==5)
	{
if(idir==0)
	    boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
  //cudaThreadSynchronize();
if(idir==1)
           boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	//boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    cudaThreadSynchronize();
	    //boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	//    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
            //boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
if(idir==2)
	    boundary2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    cudaThreadSynchronize();
	#endif
	}

	if((*p)->boundtype[field][idir][0]==0)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==1)
	boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    cudaThreadSynchronize();
	    //boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    //cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
            boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    //boundary2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    cudaThreadSynchronize();
	#endif
	}


	if((*p)->boundtype[field][idir][0]==4)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==0)
	;//boundary0_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
	;//    cudaThreadSynchronize();

if(idir==1)
	;//boundary1_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	 //   boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
      ;//      boundary2_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	   // boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	  //  cudaThreadSynchronize();
	#endif
	}

}



}



int cuboundary2(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp, struct state **d_s,  real **d_wmod,  int order,int idir,int field)
{


 dim3 dimBlock(dimblock, 1);

int numBlocks = ((dimproduct_b(*p)+numThreadsPerBlock-1)) / numThreadsPerBlock;




int i1,i2,i3;



if(((*p)->it)==-1    )
{
    //cudaMemcpy(*d_p,*p, sizeof(struct params), cudaMemcpyHostToDevice);



	if((*p)->boundtype[field][idir][0]==4)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==0)
	boundary0_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
	    cudaThreadSynchronize();

if(idir==1)
	boundary1_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	 //   boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
            boundary2_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	   // boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	  //  cudaThreadSynchronize();
	#endif
	}






	if((*p)->boundtype[field][idir][0]==3)
	{


if(idir==0)
    setboundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);
    //cudaThreadSynchronize();
if(idir==1)
    setboundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);
    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
	    setboundary2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,idir,field);
	    cudaThreadSynchronize();
	#endif
}

}
else
{

/*
	    boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
            boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);

	    cudaThreadSynchronize();
	#endif    */


	if((*p)->boundtype[field][idir][0]==5)
	{
if(idir==0)
	    boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
  //cudaThreadSynchronize();
if(idir==1)
           boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	//boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    cudaThreadSynchronize();
	    //boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	//    cudaThreadSynchronize();
	#ifdef USE_SAC_3D
            //boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
if(idir==2)
	    boundary2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    cudaThreadSynchronize();
	#endif
	}

	if((*p)->boundtype[field][idir][0]==0)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==1)
	;//boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	   ;// cudaThreadSynchronize();
	    //boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	    //cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
         ;//   boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	    //boundary2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	   ;// cudaThreadSynchronize();
	#endif
	}


	if((*p)->boundtype[field][idir][0]==4)
	{
	  //  boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
          // boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
if(idir==0)
	boundary0_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
	    cudaThreadSynchronize();

if(idir==1)
	boundary1_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	 //   boundary1_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
	 //   cudaThreadSynchronize();
	#ifdef USE_SAC_3D
if(idir==2)
            boundary2_contcd4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	   // boundary2_per_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,2,field);
	  //  cudaThreadSynchronize();
	#endif
	}

}



}









