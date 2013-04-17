#include "../include/cudapars.h"
#include "../include/iotypes.h"
#include "../include/iobparams.h"
/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "../include/smaugcukernels.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "../include/gradops_i.cuh"
#include "../include/init_user_i.cuh"


//*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd

__global__ void init_parallel(struct params *p, real *wnew, real *wmod, 
    real *dwn1, real *wd, real *wtemp, real *wtemp1, real *wtemp2)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 // int index,k;
int ni=p->n[0];
  int nj=p->n[1];
#ifdef USE_SAC_3D
  int nk=p->n[2];
#endif


// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  real *u,  *v,  *h;

   int ord;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  int i,j;
  int ip,jp;
  int ii[NDIM];
   int dimp=((p->n[0]))*((p->n[1]));

   
 #ifdef USE_SAC_3D
   int kp;
  dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
/*   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#else
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif */ 

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{
		//b[i+j*(p->n[0])]=0;

                 //Define b	

 


	//apply this special condition
	//initiate alfven wave propagtion 
	//if no initial config read

	    /*for(int f=0; f<NVAR; f++)
            { 		         
                          for(ord=0;ord<(2+3*(p->rkon==1));ord++)
                              wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=0;
	    }*/


//	 __syncthreads();

			}

        	
	 __syncthreads();


    /* #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
     
               for(int f=vel1; f<NDERV; f++)
                    wd[fencode3_i(p,ii,f)]=0.0;
     

 __syncthreads(); */



     #ifdef USE_SAC_3D
      // if((p->readini==0) && ii[0]>1 && ii[1]>1  && ii[2]>1 && ii[0]<(p->n[0])-1 && ii[1]<(p->n[1])-1 && ii[2]<(p->n[2])-1)
         if((p->readini==0) && ii[0]<(p->n[0]) && ii[1]<(p->n[1])   && ii[2]<(p->n[2])) 
     #else
      // if((p->readini==0) && ii[0]>2 && ii[1]>2 && ii[0]<(p->n[0])-3 && ii[1]<(p->n[1])-3)  //this form for OZT test???? 
     
     
     //if((p->readini==0) && ii[0]>1 && ii[1]>1  && ii[0]<(p->n[0])-1 && ii[1]<(p->n[1])-1)  //this form for OZT test???? 
        if((p->readini==0) && ii[0]<(p->n[0]) && ii[1]<(p->n[1]))  //this form for BW test  //still issue here
     #endif
	{


            #ifdef ADIABHYDRO
		    if(i> (((p->n[0])/2)-2) && i<(((p->n[0])/2)+2) && j>(((p->n[1])/2)-2) && j<(((p->n[1])/2)+2) ) 
				;//w[fencode3_i(p,ii,rho)]=1.3;
            #else
                   // init_alftest (real *w, struct params *p,int i, int j)
                   // init_alftest(w,p,i,j);
                   // init_ozttest (real *w, struct params *p,int i, int j)
                   // init_ozttest(w,p,i,j);
                   // init_bwtest(w,p,i,j);

	           //default values for positions these may be updated by the initialisation routines
                   wd[fencode3_i(p,ii,delx1)]=(p->dx[0]);
		   wd[fencode3_i(p,ii,delx2)]=(p->dx[1]);
                   wd[fencode3_i(p,ii,pos1)]=(p->xmin[0])+ii[0]*(p->dx[0]);
		   wd[fencode3_i(p,ii,pos2)]=(p->xmin[1])+ii[1]*(p->dx[1]);
                 #ifdef USE_SAC_3D
		   wd[fencode3_i(p,ii,pos3)]=(p->xmin[2])+ii[2]*(p->dx[2]);
                   wd[fencode3_i(p,ii,delx3)]=(p->dx[2]);
                 #endif

                   //init_user_i(w,p,ii);  //initilise using w field

                   //commented out because spicule problem
                   //constructed on host
                   ;//init_user_i(wmod,wd,p,ii);
           #endif

	

        }
	
	 __syncthreads();


       





     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{
        /*for(int f=energyb; f<NVAR; f++)
             if(f != rhob)
                      w[fencode3_i(p,ii,f)]=0.0;*/
        //w[fencode3_i(p,ii,b2b)]=w[fencode3_i(p,ii,b3b)];
        for(int f=rho; f<NVAR; f++)
        {               
                  //wmod[fencode3_i(p,ii,f)]=w[fencode3_i(p,ii,f)];
                  //wmod[  (((3*(1+(p->rkon)))-1)*NVAR*dimp)+fencode3_i(p,ii,f)]=w[fencode3_i(p,ii,f)];              
                  dwn1[fencode3_i(p,ii,f)]=0;

                  //initial value of ord changed to 1 ensure have correct background fields set
                  for(ord=1;ord<(2+3*(p->rkon==1));ord++)
                  {
                              //only the wmod field is used w now redundant
                              wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=wmod[fencode3_i(p,ii,f)];

                              //original version using w
                              //wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=w[fencode3_i(p,ii,f)];
                              //wmod[fencode3_i(p,ii,b2b)+ord*NVAR*dimp]=w[fencode3_i(p,ii,b3b)];
                  }
  /*        int dir=0;
 for(int field=rho;field<=rho ; field++)
if( ii[0]<4 && (p->ipe)==0  && ((p)->it)==1 && ( isnan(wmod[fencode3_i(p,ii,field)])|| wmod[fencode3_i(p,ii,field)]==0 ))
        { 
    				printf("nant %d %d %d %d %lg %lg \n",ii[0],ii[1],field,dir, wmod[fencode3_i(p,ii,rho)],wmod[fencode3_i(p,ii,field)] );
}*/
                  
        }

        for(int f=tmp1; f<NTEMP; f++)
                 wtemp[fencode3_i(p,ii,f)]=0;


	/*for(int field=rho;field<=rho ; field++)
	if(  (p->ipe)==0  && (  wmod[fencode3_i(p,ii,field)]==0 ))
		{ 
	    				printf("nanti %d %d %d %d %lg %lg \n",ii[0],ii[1],field,0, wmod[fencode3_i(p,ii,rho)],wmod[fencode3_i(p,ii,field)+dimp*NVAR] );
	}*/


}

 __syncthreads();



}

__global__ void updatemod_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, real *wtemp, real *wtemp1, real *wtemp2)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 // int index,k;
int ni=p->n[0];
  int nj=p->n[1];
#ifdef USE_SAC_3D
  int nk=p->n[2];
#endif


// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  real *u,  *v,  *h;

   int ord;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  int i,j;
  int ip,jp;
  int ii[NDIM];
   int dimp=((p->n[0]))*((p->n[1]));

   
 #ifdef USE_SAC_3D
   int kp;
  dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
/*   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#else
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif */ 

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif






 
     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{
        for(int f=rho; f<NVAR; f++)
        {               
                  for(ord=1;ord<(2+3*(p->rkon==1));ord++)
                  {
                              //wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=w[fencode3_i(p,ii,f)];
                              wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=wmod[fencode3_i(p,ii,f)];

                            // if(p->ipe==0    && f==rho)
                            //    printf("wmod,w %d %d %lg %lg\n",ii[0],ii[1],wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp],w[fencode3_i(p,ii,f)]);
 
                  }
          int dir=0;


 //for(int field=rho;field<=rho ; field++)
//if( /*ii[0]<4 &&*/ (p->ipe)==0  && /*((p)->it)==1 &&*/ (/* isnan(wmod[fencode3_i(p,ii,field)])||*/ wmod[fencode3_i(p,ii,field)]==0 ))
//        { 
//    				printf("nant %d %d %d %d %lg %lg \n",ii[0],ii[1],field,dir, wmod[fencode3_i(p,ii,rho)],wmod[fencode3_i(p,ii,field)] );
//}
                  
        }




}

 __syncthreads();



}


 //initialise grid on the gpu
 //we currently don't do this to avoid use of additional memory on GPU
//set up a temporary grid

__global__ void gridsetup_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, real *wtemp, real *wtemp1, real *wtemp2, int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 // int index,k;
int ni=p->n[0];
  int nj=p->n[1];
#ifdef USE_SAC_3D
  int nk=p->n[2];
#endif


// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  real *u,  *v,  *h;

   int ord;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  int i,j;
  int ip,jp,kp;
  int ii[NDIM];
   int dimp=((p->n[0]))*((p->n[1]));
   kp=0;
   
 #ifdef USE_SAC_3D
 
  dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
/*   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#else
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif */ 

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


     #ifdef USE_SAC_3D
       if(ii[0]>0 && ii[0]<(p->n[0]-1) && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
     {


        switch(dir)
        {

    case 0:
           wtemp2[encode3p2_i(p,ip+1,jp+1,kp+1,tmpnui)]=wd[fencode3_i(p,ii,pos1)];
    break;
    case 1:
           wtemp2[encode3p2_i(p,ip+1,jp+1,kp+1,tmpnui1)]=wd[fencode3_i(p,ii,pos2)];
    break;
    #ifdef USE_SAC_3D
           case 2:
                        wtemp2[encode3p2_i(p,ip+1,jp+1,kp+1,tmpnui2)]=wd[fencode3_i(p,ii,pos3)];
           break;
     #endif
           }
     }


        	
	 __syncthreads();




       





}




 //initialise grid on the gpu
 //we currently don't do this to avoid use of additional memory on GPU
//calculate the dx values

__global__ void setupdx_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, real *wtemp, real *wtemp1, real *wtemp2, int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 // int index,k;
int ni=p->n[0];
  int nj=p->n[1];
#ifdef USE_SAC_3D
  int nk=p->n[2];
#endif


// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  real *u,  *v,  *h;

   int ord;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  int i,j;
  int ip,jp,kp;
  int ii[NDIM];
   int dimp=((p->n[0]))*((p->n[1]));

   
 #ifdef USE_SAC_3D
 
  dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
/*   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#else
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif */ 

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

   //calculate the dx values


	    switch(dir)
	    {
		     case 0:
		     #ifdef USE_SAC_3D
		       if( ii[0]>0 && ii[0]<(p->n[0])+1 && ii[1]>0 &&  ii[1]<(p->n[1])+1 && ii[2]>0 &&  ii[2]<(p->n[2])+1)
		     #else
		       if( ii[0]>0 && ii[0]<(p->n[0])+1  && ii[1]>0 && ii[1]<(p->n[1])+1)
		     #endif
	                wd[fencode3_i(p,ii,delx1)]=0.5*(wtemp2[encode3p2_i(p,ip+1,jp,kp,tmpnui)]-wtemp2[encode3p2_i(p,ip-1,jp,kp,tmpnui)]);
		     break;
	
		     case 1:
		     #ifdef USE_SAC_3D
		       if(ii[0]>0 && ii[0]<(p->n[0])+1 && ii[1]>0 &&  ii[1]<(p->n[1])+1 && ii[2]>0 &&  ii[2]<(p->n[2])+1)
		     #else
		       if(ii[0]>0 && ii[0]<(p->n[0])+1 && ii[1]>0 && ii[1]<(p->n[1])+1)
		     #endif
			wd[fencode3_i(p,ii,delx2)]=0.5*(wtemp2[encode3p2_i(p,ip,jp+1,kp,tmpnui)]-wtemp2[encode3p2_i(p,ip,jp-1,kp,tmpnui)]);
		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:

		       if(ii[0]>0 && ii[0]<(p->n[0])+1 && ii[1]>0 && ii[1]<(p->n[1])+1 && ii[2]>0 && ii[2]<(p->n[2])+1)
			wd[fencode3_i(p,ii,delx3)]=0.5*(wtemp2[encode3p2_i(p,ip,jp,kp+1,tmpnui)]-wtemp2[encode3p2_i(p,ip,jp,kp-1,tmpnui)]);
		     break;			
		     #endif
	     }
     
        	
	 __syncthreads();







       





}

 //initialise grid on the gpu
 //we currently don't do this to avoid use of additional memory on GPU
//intialise temporrary matrix needs t be completed
__global__ void zerotempv_parallel(struct params *p, real *w, real *wnew, real *wmod, 
real *dwn1,  real *wd, real *wtemp, real *wtemp1, real *wtemp2,  int dir)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  const int blockdim=blockDim.x;
  const int SZWT=blockdim;
  const int SZWM=blockdim*NVAR;
  int tid=threadIdx.x;
  real maxt=0,max3=0, max1=0;
  int i,j,iv;
  int is,js;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];


  
   int ip,jp;



  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif  

int bfac1,bfac2,bfac3;
//int bfac1=(field==rho || field>mom2)+(field>rho && field<energy);
//int bfac2= (field==rho || field>mom2);
//int bfac3=(field>rho && field<energy);
//int shift=order*NVAR*dimp;




//init temp1 and temp2 to zero 
//the compute element initialising n[0] or n[1] element must do +1 and +2
//this is because we fit the problem geometrically to nixnj elements 

     ii[0]=ip;
     ii[1]=jp;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
           k=ii[2];
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
    //set viscosities
   //if(i<((p->n[0])) && j<((p->n[1])))
   {


        for(int f=d1; f<=d3; f++)
     #ifdef USE_SAC_3D
                 wtemp2[encode3p2_i(p,ii[0],ii[1],ii[2],tmpnui)]=0;
     #else
                 wtemp2[encode3p2_i(p,ii[0],ii[1],k,tmpnui)]=0;
     #endif

      if(i==((p->n[0])-1))
      {
        wtemp2[encode3p2_i(p,i+1,j,k,tmpnui)]=0;
        wtemp2[encode3p2_i(p,i+2,j,k,tmpnui)]=0;
      }
      if(j==((p->n[1])-1))
      {
          wtemp2[encode3p2_i(p,i,j+1,k,tmpnui)]=0;
          wtemp2[encode3p2_i(p,i,j+2,k,tmpnui)]=0;
      }

     #ifdef USE_SAC_3D
      if(k==((p->n[2])-1))
      {
          wtemp2[encode3p2_i(p,i,j,k+1,tmpnui)]=0;
          wtemp2[encode3p2_i(p,i,j,k+2,tmpnui)]=0;
      }

     #endif
      if(j==((p->n[1])-1)  && i==((p->n[0])-1))
      {
          for(int di=0; di<2; di++)
             for(int dj=0; dj<2; dj++)
                   wtemp2[encode3p2_i(p,i+1+di,j+1+dj,k,tmpnui)]=0;
      }
     #ifdef USE_SAC_3D
      if(i==((p->n[0])-1)  && k==((p->n[2])-1))
      {
          for(int di=0; di<2; di++)
             for(int dk=0; dk<2; dk++)
                   wtemp2[encode3p2_i(p,i+1+di,j,k+1+dk,tmpnui)]=0;
      }
      #endif

    

     #ifdef USE_SAC_3D
      if(j==((p->n[1])-1)  && k==((p->n[2])-1))
      {
          for(int dk=0; dk<2; dk++)
             for(int dj=0; dj<2; dj++)
                   wtemp2[encode3p2_i(p,i,j+1+dj,k+1+dk,tmpnui)]=0;
      }
      #endif

     #ifdef USE_SAC_3D
      if(i==((p->n[0])-1) && j==((p->n[1])-1)  && k==((p->n[2])-1))
      {
          for(int dk=0; dk<2; dk++)
             for(int dj=0; dj<2; dj++)
               for(int di=0; di<2; di++)
                   wtemp2[encode3p2_i(p,i+1+di,j+1+dj,k+1+dk,tmpnui)]=0;
      }
      #endif

   }

}



__device__ __host__
int encodempiw (struct params *p,int ix, int iy, int iz, int field,int bound,int dim) {
  #ifdef USE_SAC_3D
    return (dim*(    4*NVAR*(         ((p->n[0])*(p->n[1]))+((p->n[1])*(p->n[2]))+((p->n[0])*(p->n[2]))   )           )+4*field*(         ((p->n[0])*(p->n[1]))+((p->n[1])*(p->n[2]))+((p->n[0])*(p->n[2]))   )+
bound*(         (dim==2)*((p->n[0])*(p->n[1]))   +  (dim==0)*((p->n[1])*(p->n[2]))  +   (dim==1)*((p->n[0])*(p->n[2]))    )+   (  (ix+iz*(p->n[0]))*(dim==1)+(iy+iz*(p->n[1]))*(dim==0)+(iz+ix*(p->n[2]))*(dim==2)    ));
  #else
    return (dim*(4*NVAR*((p->n[0])+(p->n[1])))+4*field*((p->n[0])+(p->n[1]))+bound*((dim==1)*(p->n[0])+(dim==0)*(p->n[1]))  +   (ix*(dim==1)+iy*(dim==0)));
  #endif
}

__device__ __host__
int encodempiw0 (struct params *p,int ix, int iy, int iz, int field,int bound) {
  #ifdef USE_SAC_3D
    return (4*field*(         ((p->n[1])*(p->n[2]))   )+
bound*(            +  ((p->n[1])*(p->n[2]))      )+   (  (iy+iz*(p->n[1]))    ));
  #else
    return (   4*field*(p->n[1]) +bound*((p->n[1]))  +   (iy)   );
  #endif
}


__device__ __host__
int encodempiw1 (struct params *p,int ix, int iy, int iz, int field,int bound) {
  #ifdef USE_SAC_3D
    return (4*field*(         ((p->n[0])*(p->n[2]))   )+
bound*(            +  ((p->n[0])*(p->n[2]))      )+   (  (ix+iz*(p->n[0]))    ));
  #else
    return (4*field*(p->n[0]) +bound*((p->n[0]))  +   (ix));
  #endif
}

__device__ __host__
int encodempiw2 (struct params *p,int ix, int iy, int iz, int field,int bound) {
  #ifdef USE_SAC_3D
    return (4*field*(         ((p->n[0])*(p->n[1]))   )+
bound*(            +  ((p->n[0])*(p->n[1]))      )+   (  (ix+iy*(p->n[0]))    ));
  #endif
}


__device__ __host__
int encodempivisc (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (dim*(    2*(         (((p->n[0])+2)*((p->n[1])+2))+(((p->n[1])+2)*((p->n[2])+2))+(((p->n[0])+2)*((p->n[2])+2))   )           )+
bound*(         (dim==2)*(((p->n[0])+2)*((p->n[1])+2))   +  (dim==0)*(((p->n[1])+2)*((p->n[2])+2))  +   (dim==1)*(((p->n[0])+2)*((p->n[2])+2))    )+   (  (ix+iz*((p->n[0])+2))*(dim==1)+(iy+iz*((p->n[1])+2))*(dim==0)+(iz+ix*((p->n[2])+2))*(dim==2)    ));
  #else
    return (   dim*(2*(  ((p->n[0])+2)+((p->n[1])+2)   ))      +bound*(    (dim==1)*((p->n[0])+2)+(dim==0)*((p->n[1])+2)  )  +   (ix*(dim==1)+iy*(dim==0))     );
  #endif
}


__device__ __host__
int encodempivisc0 (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (
bound*(           (((p->n[1])+2)*((p->n[2])+2))      )+   (  (iy+iz*((p->n[1])+2))    ));
  #else
    return (   bound*(    ((p->n[1])+2)  )  +   iy     );
  #endif
}


__device__ __host__
int encodempivisc1 (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (
bound*(           (((p->n[0])+2)*((p->n[2])+2))      )+   (  (ix+iz*((p->n[0])+2))    ));
  #else
    return (   bound*(    ((p->n[0])+2)  )  +   ix     );
  #endif
}

__device__ __host__
int encodempivisc2 (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (
bound*(           (((p->n[0])+2)*((p->n[1])+2))      )+   (  (ix+iy*((p->n[0])+2))    ));
  #endif
}


#ifdef USE_MPI

     __device__ __host__ void mpiwtogpu(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int *ii, int var, int dim)
    {

             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;




                if((i==0 || i==1) && dim==0)
                {              
                    bound=i;
                    d_w[encode3_i(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                   // d_wmod[encode3_i(p,i,j,k,var)]=d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)];
                    // if(var==4 && ((p)->ipe)==0)                        
                    //    printf(" %d %d %d %d actual %d  mpi data%d %g\n",i,j,bound,dim,var,encodempiw0(p,i,j,k,var,bound),d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);   


 

     
      
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_w[encode3_i(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                  //  d_wmod[encode3_i(p,i,j,k,var)]=d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)];    


                          // if(p->ipe==0    && var==rho && dim==0 )
                          //   {
                             // for(int bound=0;bound<=1;bound++)
                             //   printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiw0[encodempiw0(p,i,j,k,var,bound)]);
			     	//printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                           //  }
          
                }

              

                if((j==0 || j==1) && dim==1)              
                {              
                    bound=j;
                    d_w[encode3_i(p,i,j,k,var)]=d_mpiw1[encodempiw1(p,i,j,k,var,bound)];
                   // d_wmod[encode3_i(p,i,j,k,var)]=d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)];              
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1)               
                {
                   bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_w[encode3_i(p,i,j,k,var)]=d_mpiw1[encodempiw1(p,i,j,k,var,bound)];
                  //  d_wmod[encode3_i(p,i,j,k,var)]=d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)];              
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                    d_w[encode3_i(p,i,j,k,var)]=d_mpiw2[encodempiw2(p,i,j,k,var,bound)];
                  //  d_wmod[encode3_i(p,i,j,k,var)]=d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)];              
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_w[encode3_i(p,i,j,k,var)]=d_mpiw2[encodempiw2(p,i,j,k,var,bound)];
                  //  d_wmod[encode3_i(p,i,j,k,var)]=d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)];              
                }

     #endif
 
// if( (p->ipe)==0  && ((p)->it)==0 && (isnan(d_wmod[fencode3_i(p,ii,rho)]) || d_wmod[fencode3_i(p,ii,rho)]==0      ))
//        { 
//    	printf("nant %d %d %d %lg\n",ii[0],ii[1], dim, d_wmod[fencode3_i(p,ii,rho)] );
//}

    }

     __device__ __host__ void mpiwmodtogpu(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int *ii, int var, int dim, int order)
    {

             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;

 #ifdef USE_SAC_3D
	int ntot=((p->n[0]))*((p->n[1]))*((p->n[2]))*NVAR;
 #else
	int ntot=((p->n[0]))*((p->n[1]))*NVAR;
 #endif
                //remember only update the boundaries if they are mpiupper boundaries 
                //or an mpi period 


                if((i==0 || i==1) && dim==0  && ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))    )
                {              
                    bound=i;
                    d_wmod[order*ntot+encode3_i(p,i,j,k,var)]=d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)];
                    // if(var==0 && ((p)->ipe)==0)                        
                    //    printf(" %d %d %d %d actual mpi data %d %g\n",i,j,bound,dim,encodempiw0(p,i,j,k,var,bound),d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]); 

       
      
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0  && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2)) )               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_wmod[order*ntot+encode3_i(p,i,j,k,var)]=d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)];  

                             //if(p->ipe==0    && var==rho && dim==0 )
                             //{
                             // for(int bound=0;bound<=1;bound++)
                             //   printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);
			     //	printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                             //}
          
                }

              

                if((j==0 || j==1) && dim==1  && ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))   )              
                {              
                    bound=j;
                    d_wmod[order*ntot+encode3_i(p,i,j,k,var)]=d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]; 



             
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1   && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2))    )               
                {
                   bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_wmod[order*ntot+encode3_i(p,i,j,k,var)]=d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]; 

  

             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2   && ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))   )              
                {              
                    bound=k;
                    d_wmod[order*ntot+encode3_i(p,i,j,k,var)]=d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)];              
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2   && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2))   )               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_wmod[order*ntot+encode3_i(p,i,j,k,var)]=d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)];              
                }

     #endif
 
// if( (p->ipe)==0  && ((p)->it)==0 && (isnan(d_wmod[fencode3_i(p,ii,rho)]) || d_wmod[fencode3_i(p,ii,rho)]==0      ))
//        { 
//    	printf("nant %d %d %d %lg\n",ii[0],ii[1], dim, d_wmod[fencode3_i(p,ii,rho)] );
//}

    }



     __device__ __host__ void mpiwdtogpu(struct params *p,real *d_wd,real *d_mpiw0,real *d_mpiw1,real *d_mpiw2,int *ii, int var, int dim)
    {

             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;

                //remember only update the boundaries if they are mpiupper boundaries 
                //or an mpi period 



 
                if((i==0 || i==1) && dim==0   &&  ((p->mpilowerb[dim])==1))
                {              
                    bound=i;
                    d_wd[encode3_i(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                    // if(var==4 && ((p)->ipe)==0)                        
                    //    printf(" %d %d %d %d actual %d  mpi data%d %g\n",i,j,bound,dim,var,encodempiw0(p,i,j,k,var,bound),d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);        
      
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0  &&  ((p->mpiupperb[dim])==1))               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_wd[encode3_i(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                }

              

                if((j==0 || j==1) && dim==1   &&  ((p->mpilowerb[dim])==1))              
                {              
                    bound=j;
                    d_wd[encode3_i(p,i,j,k,var)]=d_mpiw1[encodempiw1(p,i,j,k,var,bound)];
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1  &&  ((p->mpiupperb[dim])==1))               
                {
                   bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_wd[encode3_i(p,i,j,k,var)]=d_mpiw1[encodempiw1(p,i,j,k,var,bound)];
                }


               /* if((i==0 || i==1) && (j==0 || j==1))
                {              
                    bound=i;
                    d_wd[encode3_i(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                    // if(var==4 && ((p)->ipe)==0)                        
                    //    printf(" %d %d %d %d actual %d  mpi data%d %g\n",i,j,bound,dim,var,encodempiw0(p,i,j,k,var,bound),d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);        
      
                }
                if((( j>=((p->n[1])-2)   ))  && (( i>=((p->n[0])-2)   )))               
                {
                   bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_wd[encode3_i(p,i,j,k,var)]=d_mpiw1[encodempiw0(p,i,j,k,var,bound)];
                }*/
                



       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                    d_wd[encode3_i(p,i,j,k,var)]=d_mpiw2[encodempiw2(p,i,j,k,var,bound)];
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_wd[encode3_i(p,i,j,k,var)]=d_mpiw2[encodempiw2(p,i,j,k,var,bound)];
                }

     #endif
 


    }




__device__ __host__ void   mpivisctogpu(struct params *p,real *d_wtemp2,real *d_gmpivisc0,real *d_gmpivisc1,real *d_gmpivisc2,int *ii,  int dim)
{
                                
               int i,j,k,bound,var;
              var=0;
i=ii[0];
j=ii[1];
k=0;
                //remember only update the boundaries if they are mpiupper boundaries 
                //or an mpi period 
 
                if((i==0 ) && dim==0 /* && ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)
                {              
                    bound=i;
                    d_wtemp2[encode3p2_i(p,i,j,k,var)]=d_gmpivisc0[encodempivisc0(p,i,j,k,bound,dim)];
                    
                }
                else if((( i==((p->n[0])+1)   ))  && dim==0  /* && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2)) */ )               
                {
                    bound=1;
                    d_wtemp2[encode3p2_i(p,i,j,k,var)]=d_gmpivisc0[encodempivisc0(p,i,j,k,bound,dim)];
                }

              

                if((j==0) && dim==1 /* && ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)              
                {              
                    bound=j;
                    d_wtemp2[encode3p2_i(p,i,j,k,var)]=d_gmpivisc1[encodempivisc1(p,i,j,k,bound,dim)];
                }            
                 else if((( j==((p->n[1])+1)   ))  && dim==1   /*&& ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)               
                {
                    bound=1;
                    d_wtemp2[encode3p2_i(p,i,j,k,var)]=d_gmpivisc1[encodempivisc1(p,i,j,k,bound,dim)];
             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 ) && dim==2  /*&& ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)              
                {              
                    bound=k;
                    d_wtemp2[encode3p2_i(p,i,j,k,var)]=d_gmpivisc2[encodempivisc2(p,i,j,k,bound,dim)];
                }        
                 else if(((k==((p->n[2])+1)   ))  && dim==2   /* && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)               
                {
                    bound=1;
                    d_wtemp2[encode3p2_i(p,i,j,k,var)]=d_gmpivisc2[encodempivisc2(p,i,j,k,bound,dim)];
                }

     #endif
                               
                                
}

__device__ __host__ void   gputompivisc(struct params *p,real *d_wtemp2,real *d_gmpivisc0,real *d_gmpivisc1,real *d_gmpivisc2,int *ii,  int dim)
{
                                
              int i,j,k,bound,var;
              var=0;
i=ii[0];
j=ii[1];
k=0;
 
 
                if((i==0 ) && dim==0)
                {              
                    bound=i;
                    d_gmpivisc0[encodempivisc0(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_i(p,i,j,k,var)];
                    
                }
                else if((( i==((p->n[0])+1)   ))  && dim==0)               
                {
                    bound=1;
                    d_gmpivisc0[encodempivisc0(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_i(p,i,j,k,var)];
                }

              

                if((j==0) && dim==1)              
                {              
                    bound=j;
                    d_gmpivisc1[encodempivisc1(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_i(p,i,j,k,var)];
                }            
                 else if((( j==((p->n[1])+1)   ))  && dim==1)               
                {
                    bound=1;
                    d_gmpivisc1[encodempivisc1(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_i(p,i,j,k,var)];
             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 ) && dim==2)              
                {              
                    bound=k;
                    d_gmpivisc2[encodempivisc2(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_i(p,i,j,k,var)];
                }        
                 else if(((k==((p->n[2])+1)   ))  && dim==2)               
                {
                    bound=1;
                    d_gmpivisc2[encodempivisc2(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_i(p,i,j,k,var)];
                }

     #endif
                               
                                
}

     __device__ __host__ void gputompiw(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int *ii, int var, int dim)
    {
             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;

/*  for(int field=rho;field<=rho ; field++)
if(  (p->ipe)==0  && ((p)->it)==1 && ( isnan(d_wmod[fencode3_i(p,ii,field)])|| d_wmod[fencode3_i(p,ii,field)]==0 ))
        { 
    				printf("nant %d %d %d %d %lg %lg \n",ii[0],ii[1],field,dim,d_wmod[fencode3_i(p,ii,rho)],d_wmod[fencode3_i(p,ii,field)] );
}*/
 
                if((i==0 || i==1) && dim==0)
                {              
                    bound=i;
                    d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_w[encode3_i(p,i+2,j,k,var)];
                   // d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[encode3_i(p,i+2,j,k,var)];



              
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_w[encode3_i(p,i-2,j,k,var)];
                   // d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[encode3_i(p,i-2,j,k,var)];  




                }

              

                if((j==0 || j==1) && dim==1)              
                {              
                    bound=j;
                    d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_w[encode3_i(p,i,j+2,k,var)];
                   // d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]=d_wmod[encode3_i(p,i,j+2,k,var)];   








			/*  if( var==0 && (p)->ipe==3 && ((p)->it)==2 && bound==3)
			    {
				 //printf("ipe3 mpiwmod \n");
				 //for(int iii=0; iii<4*((p)->n[0]);iii++)
				     printf(" %lg %d \n",d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)],encodempiw1(p,i,j,k,var,bound));
				 //printf("\n");
			     }*/


           
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1)               
                {
                    bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_w[encode3_i(p,i,j-2,k,var)];
                   // d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]=d_wmod[encode3_i(p,i,j-2,k,var)];  

                          /*   if(p->ipe==0    && var==rho && dim==1 )
                             {
                             // for(int bound=0;bound<=1;bound++)
                                printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiw0[encodempiw0(p,i,j,k,var,bound)]);
			     	//printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                             }*/



             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                    d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_w[encode3_i(p,i,j,k+2,var)];
                   // d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)]=d_wmod[encode3_i(p,i,j,k+2,var)];              
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_w[encode3_i(p,i,j,k-2,var)];
                   // d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)]=d_wmod[encode3_i(p,i,j,k-2,var)];               
                }

     #endif



    /*if((p)->ipe==3 && ((p)->it)==2  && i==((p->n[0])-1) && j==((p->n[1])-1))
    {
         printf("ipe3 mpiwmod \n");
         for(int iii=0; iii<4*((p)->n[0]);iii++)
             printf(" %lg ",d_mpiwmod1[iii]);
         printf("\n");
     }*/

 
 
 }


     __device__ __host__ void gputompiwmod(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int *ii, int var, int dim, int order)
    {
             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;

 #ifdef USE_SAC_3D
	int ntot=((p->n[0]))*((p->n[1]))*((p->n[2]))*NVAR;
 #else
	int ntot=((p->n[0]))*((p->n[1]))*NVAR;
 #endif

/*  for(int field=rho;field<=rho ; field++)
if(  (p->ipe)==0  && ((p)->it)==1 && ( isnan(d_wmod[fencode3_i(p,ii,field)])|| d_wmod[fencode3_i(p,ii,field)]==0 ))
        { 
    				printf("nant %d %d %d %d %lg %lg \n",ii[0],ii[1],field,dim,d_wmod[fencode3_i(p,ii,rho)],d_wmod[fencode3_i(p,ii,field)] );
}*/
 
                if((i==0 || i==1) && dim==0)
                {              
                    bound=i;
                    d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_i(p,i+2,j,k,var)];
                   // d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[encode3_i(p,i+2,j,k,var)];

                            // if(p->ipe==0    && var==rho && dim==0 )
                             //{
                             // for(int bound=0;bound<=1;bound++)
                              //  printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);
			     	//printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                             //}
              
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    //d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_w[encode3_i(p,i-2,j,k,var)];
                    d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_i(p,i-2,j,k,var)];  

                             //if(p->ipe==0    && var==rho && dim==0 )
                             //{
                             // for(int bound=0;bound<=1;bound++)
                             //   printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);
			     	//printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                             //}





                }

              

                if((j==0 || j==1) && dim==1)              
                {              
                    bound=j;
                    //d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_w[encode3_i(p,i,j+2,k,var)];
                    d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_i(p,i,j+2,k,var)];   






			/*  if( var==0 && (p)->ipe==3 && ((p)->it)==2 && bound==3)
			    {
				 //printf("ipe3 mpiwmod \n");
				 //for(int iii=0; iii<4*((p)->n[0]);iii++)
				     printf(" %lg %d \n",d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)],encodempiw1(p,i,j,k,var,bound));
				 //printf("\n");
			     }*/


           
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1)               
                {
                    bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                   // d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_w[encode3_i(p,i,j-2,k,var)];
                    d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_i(p,i,j-2,k,var)];  





             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                   // d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_w[encode3_i(p,i,j,k+2,var)];
                    d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_i(p,i,j,k+2,var)];              
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    //d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_w[encode3_i(p,i,j,k-2,var)];
                    d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_i(p,i,j,k-2,var)];               
                }

     #endif



    /*if((p)->ipe==3 && ((p)->it)==2  && i==((p->n[0])-1) && j==((p->n[1])-1))
    {
         printf("ipe3 mpiwmod \n");
         for(int iii=0; iii<4*((p)->n[0]);iii++)
             printf(" %lg ",d_mpiwmod1[iii]);
         printf("\n");
     }*/

 
 
 }





     __device__ __host__ void gputompiwd(struct params *p,real *d_wd,real *d_mpiw0,real *d_mpiw1,real *d_mpiw2,int *ii, int var, int dim)
    {
             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;
 
 
                if((i==0 || i==1) && dim==0)
                {              
                    bound=i;
                    d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_wd[encode3_i(p,i+2,j,k,var)];
                 }
                else if((( i>=((p->n[0])-2)   ))  && dim==0)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_wd[encode3_i(p,i-2,j,k,var)];
                  }

              

                if((j==0 || j==1) && dim==1)              
                {              
                    bound=j;
                    d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_wd[encode3_i(p,i,j+2,k,var)];
 

			/*  if( var==0 && (p)->ipe==3 && ((p)->it)==2 && bound==3)
			    {
				 //printf("ipe3 mpiwmod \n");
				 //for(int iii=0; iii<4*((p)->n[0]);iii++)
				     printf(" %lg %d \n",d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)],encodempiw1(p,i,j,k,var,bound));
				 //printf("\n");
			     }*/


           
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1)               
                {
                    bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_wd[encode3_i(p,i,j-2,k,var)];
                 }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                    d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_wd[encode3_i(p,i,j,k+2,var)];
                 }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_wd[encode3_i(p,i,j,k-2,var)];
                 }

     #endif



    /*if((p)->ipe==3 && ((p)->it)==2  && i==((p->n[0])-1) && j==((p->n[1])-1))
    {
         printf("ipe3 mpiwmod \n");
         for(int iii=0; iii<4*((p)->n[0]);iii++)
             printf(" %lg ",d_mpiwmod1[iii]);
         printf("\n");
     }*/

 
 
 }






__global__ void  mpiwtogpu_parallel(struct params *p,real *d_w, real *d_wmod, real *d_mpiw0, real *d_mpiwmod0, real *d_mpiw1, real *d_mpiwmod1, real *d_mpiw2, real *d_mpiwmod2, int idir)
{

int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

#ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


//int shift=order*NVAR*dimp;

//if(iindex==0)
//         printf("in mpiwtogpu\n");

     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
            for( f=rho; f<NVAR; f++)
     #else
     //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	  for( f=rho; f<NVAR; f++)
     #endif     
         #ifdef USE_SAC_3D
           if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
         #else
           if(i<((p->n[0])) && j<((p->n[1])))
         #endif     
{      
    // if(iindex==0)
    //     printf("calling  mpiwtogpu %d %d\n",dim,f);

                 mpiwtogpu(p,d_w,d_wmod,d_mpiw0,d_mpiwmod0,d_mpiw1,d_mpiwmod1,d_mpiw2,d_mpiwmod2,iia,f,idir);

}


 __syncthreads();

           
               
}



__global__ void  mpiwmodtogpu_parallel(struct params *p,real *d_w, real *d_wmod, real *d_mpiw0, real *d_mpiwmod0, real *d_mpiw1, real *d_mpiwmod1, real *d_mpiw2, real *d_mpiwmod2, int idir, int order)
{

int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

//int var,bound;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

#ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


//int shift=order*NVAR*dimp;

//if(iindex==0)
//         printf("in mpiwtogpu\n");

     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
            for( f=rho; f<NVAR; f++)
     #else
     //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	  for( f=rho; f<NVAR; f++)
     #endif     
         #ifdef USE_SAC_3D
           if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
         #else
           if(i<((p->n[0])) && j<((p->n[1])))
         #endif     
{      
    // if(iindex==0)
    //     printf("calling  mpiwtogpu %d %d\n",dim,f);

                    //bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    /*bound=i;
                    var=f;
                    d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[encode3_i(p,i+2,j,k,var)];


                 if( f==rho && (p->ipe)==0  && (  d_wmod[fencode3_i(p,iia,f)]==0 )   && i<=1)
                       if(idir==0)
    				printf("nani0 %d %d  %lg  \n",iia[0],iia[1], d_mpiwmod0[encodempiw0(p,i,j,k,f,bound)] );*/



                 mpiwmodtogpu(p,d_w,d_wmod,d_mpiw0,d_mpiwmod0,d_mpiw1,d_mpiwmod1,d_mpiw2,d_mpiwmod2,iia,f,idir,order);



                // if( f==rho && (p->ipe)==0  && (  d_wmod[fencode3_i(p,iia,f)]==0 ))
                //       if(idir==0)
    		//		printf("nani0 %d %d  %lg %lg \n",iia[0],iia[1], d_wmod[fencode3_i(p,iia,rho)],d_wmod[fencode3_i(p,iia,f)+dimp*NVAR] );

}


 __syncthreads();

           
               
}




__global__ void  mpiwdtogpu_parallel(struct params *p,  int dir, int var, real *d_wd, real *d_mpiw0,  real *d_mpiw1,  real *d_mpiw2, int idir)
{

int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

#ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


//int shift=order*NVAR*dimp;

//if(iindex==0)
//         printf("in mpiwtogpu\n");

     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
            for( f=pos1; f<=delx3; f++)
     #else
     //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	  for( f=pos1; f<=delx2; f++)
     #endif     
         #ifdef USE_SAC_3D
           if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
         #else
           if(i<((p->n[0])) && j<((p->n[1])))
         #endif     
{      
    // if(iindex==0)
    //     printf("calling  mpiwtogpu %d %d\n",dim,f);

                 mpiwdtogpu(p,d_wd,d_mpiw0,d_mpiw1,d_mpiw2,iia,f,idir);

}


 __syncthreads();

           
               
}



     __global__ void gputompiwmod_parallel(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int order, int idir)
    {

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


//int shift=order*NVAR*dimp;


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
           for( f=rho; f<NVAR; f++)
     #else
           //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	   for( f=rho; f<NVAR; f++)
     #endif
             {
            
         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{

 

                  gputompiwmod(p,d_w,d_wmod,d_mpiw0,d_mpiwmod0,d_mpiw1,d_mpiwmod1,d_mpiw2,d_mpiwmod2,iia,f,idir,order);

                            // if(p->ipe==0    && f==rho && idir==0 )
                            //  for(int bound=0;bound<=1;bound++)
                            //    printf("mpiw0 %d %d %d %d %lg %lg\n",idir,bound,iia[0],iia[1],d_mpiwmod0[encodempiw0(p,i,j,k,f,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,f,bound)]);

	}

               }





 __syncthreads();

}

     __global__ void gputompiw_parallel(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int order, int idir)
    {

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


//int shift=order*NVAR*dimp;


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
           for( f=rho; f<NVAR; f++)
     #else
           //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	   for( f=rho; f<NVAR; f++)
     #endif
             {
            
         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{

 

                  gputompiw(p,d_w,d_wmod,d_mpiw0,d_mpiwmod0,d_mpiw1,d_mpiwmod1,d_mpiw2,d_mpiwmod2,iia,f,idir);

                             //if(p->ipe==0    && f==rho && idir==0 )
                             // for(int bound=0;bound<=1;bound++)
                             //   printf("mpiw0 %d %d %d %d %lg %lg\n",idir,bound,iia[0],iia[1],d_mpiw0[encodempiw0(p,i,j,k,f,bound)],d_mpiw1[encodempiw1(p,i,j,k,f,bound)]);

	}

               }





 __syncthreads();

}


     __global__ void gputompiwd_parallel(struct params *p,real *d_wd,real *d_mpiw0,real *d_mpiw1,real *d_mpiw2,int order, int idir)
    {

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


//int shift=order*NVAR*dimp;


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
           for( f=pos1; f<=delx3; f++)
     #else
        //   for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	   for( f=pos1; f<=delx2; f++)
     #endif
             {
            
         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{

 

                  gputompiwd(p,d_wd,d_mpiw0,d_mpiw1,d_mpiw2,iia,f,idir);


	}

               }





 __syncthreads();

}


     __global__ void gputompivisc_parallel(struct params *p,real *d_wtemp2,real *d_gmpivisc0,real *d_gmpivisc1,real *d_gmpivisc2)
     {
               
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/((nj+2)*(ni+2));
   jp=(iindex-(kp*((nj+2)*(ni+2))))/(ni+2);
   ip=iindex-(kp*(nj+2)*(ni+2))-(jp*(ni+2));
#else
    jp=iindex/(ni+2);
   ip=iindex-(jp*(ni+2));
#endif     


//int shift=order*NVAR*dimp;


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];

     #else

     #endif
           for(dim=0; dim<NDIM;dim++)
             {
            
         #ifdef USE_SAC_3D
      if(i<(((p->n[0])+2)) && j<(((p->n[1])+2))  && k<(((p->n[2])+2)))
     #else
       if(i<(((p->n[0])+2)) && j<(((p->n[1])+2)))
     #endif           
	{

 

                  gputompivisc(p,d_wtemp2,d_gmpivisc0,d_gmpivisc1,d_gmpivisc2,iia,dim);

	}

               }

 __syncthreads();
              
               }    
     
     
    __global__ void  mpivisctogpu_parallel(struct params *p,real *d_wtemp2,real *d_gmpivisc0,real *d_gmpivisc1,real *d_gmpivisc2)
    {
               
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/((nj+2)*(ni+2));
   jp=(iindex-(kp*((nj+2)*(ni+2))))/(ni+2);
   ip=iindex-(kp*(nj+2)*(ni+2))-(jp*(ni+2));
#else
    jp=iindex/(ni+2);
   ip=iindex-(jp*(ni+2));
#endif     


//int shift=order*NVAR*dimp;


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];

     #else

     #endif
           for(dim=0; dim<NDIM;dim++)
             {
            
         #ifdef USE_SAC_3D
      if(i<(((p->n[0])+2)) && j<(((p->n[1])+2))  && k<(((p->n[2])+2)))
     #else
       if(i<(((p->n[0])+2)) && j<(((p->n[1])+2)))
     #endif           
	{

 

                  mpivisctogpu(p,d_wtemp2,d_gmpivisc0,d_gmpivisc1,d_gmpivisc2,iia,dim);

	}

               }

 __syncthreads();
               
               
}

#endif



/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_i(char *label)
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

int cusync(struct params **p)
{

  #ifdef USE_GPUD
     
         for(int igid=0; igid<((*p)->npe); igid++)
         {
                (*p)->ipe=igid;
                cudaSetDevice((*p)->gpid[igid]) ;
                
  #endif
  cudaThreadSynchronize();
  #ifdef USE_GPUD
                 (*p)->ipe=0;
                 cudaSetDevice((*p)->gpid[0]) ;
          }
  #endif
  return 0;
}

int cusetgpu(struct params **p)
{
  #ifdef USE_GPUD
    if(((*p)->ipe)==-1)
    {
         for(int igid=0; igid<((*p)->npe); igid++)
                (*p)->gpid[igid]=igid ;
    }
    else
      cudaSetDevice((*p)->gpid[(*p)->ipe]) ;
                
  #endif
 
  return 0;
}

int cuinit(struct params **p, struct bparams **bp, real **wmod,real **wnew, real **wd, struct state **state, struct params **d_p, struct bparams **d_bp, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{



/////////////////////////////////////
  // (1) initialisations:
  //     - perform basic sanity checks
  //     - set device
  /////////////////////////////////////
  int deviceCount;
  int dir;
 /* cudaGetDeviceCount(&deviceCount);
   
 // if (deviceCount == 0)
 // {
 //   fprintf(stderr, "Sorry, no CUDA device fount");
 //   return 1;
//  }

  #ifdef USE_MPI
     int lipe=(*p)->ipe;
     int gpugid=lipe/4;
     selectedDevice=lipe-4*gpugid;
  #endif
  if (selectedDevice >= deviceCount)
  {
    fprintf(stderr, "Choose device ID between 0 and %d\n", deviceCount-1);
    return 1;
  }


        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, selectedDevice);
        if (deviceProp.major < 1) {
            fprintf(stderr, "gpuDeviceInit(): GPU device does not support CUDA.\n");
            exit(-1);                                                  \
        }

        cudaSetDevice(selectedDevice) ;
        printf("> gpuDeviceInit() CUDA device [%d]: %s %s\n", selectedDevice, deviceProp.name, getenv("HOSTNAME"));



  cudaSetDevice(selectedDevice);
  printf("device count %d selected %d\n", deviceCount,selectedDevice);
  checkErrors_i("initialisations");*/
  
	// Build empty u, v, b matrices

  printf("in cuinit\n");
 // real *adb;
  real *adw, *adwnew;
  struct params *adp;
  struct bparams *adbp;
  struct state *ads;


 
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif  
if((*p)->mode != 3)
{
	if(((*p)->rkon)==1)
	  cudaMalloc((void**)d_wmod, 6*NVAR*dimp*sizeof(real));
	else
	  cudaMalloc((void**)d_wmod, 3*NVAR*dimp*sizeof(real));

	  cudaMalloc((void**)d_dwn1, NVAR*dimp*sizeof(real));
	  cudaMalloc((void**)d_wd, NDERV*dimp*sizeof(real));
	  cudaMalloc((void**)d_wtemp, NTEMP*dimp*sizeof(real));


	  #ifdef USE_SAC
	  cudaMalloc((void**)d_wtemp1, NTEMP1*(((*p)->n[0])+1)* (((*p)->n[1])+1)*sizeof(real));
	  cudaMalloc((void**)d_wtemp2, NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)*sizeof(real));
	  #endif
	  #ifdef USE_SAC_3D
	  cudaMalloc((void**)d_wtemp1, NTEMP1*(((*p)->n[0])+1)* (((*p)->n[1])+1)* (((*p)->n[2])+1)*sizeof(real));
	  cudaMalloc((void**)d_wtemp2, NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2)*sizeof(real));
	  #endif

	  //cudaMalloc((void**)&adw, NVAR*dimp*sizeof(real));
	  //cudaMalloc((void**)&adwnew, NVAR*dimp*sizeof(real));

	  cudaMalloc((void**)&adbp, sizeof(struct bparams));
	  cudaMalloc((void**)&adp, sizeof(struct params));
	  cudaMalloc((void**)&ads, sizeof(struct state));
	 // checkErrors_i("memory allocation");

	printf("ni is %d\n",(*p)->n[1]);

	   // *d_b=adb;
	    *d_bp=adbp;
	    *d_p=adp;
	    //*d_w=adw;
	    //*d_wnew=adwnew;
	    *d_state=ads;

	     
	//printf("allocating %d %d %d %d\n",dimp,(*p)->n[0],(*p)->n[1],(*p)->n[2]);
	printf("allocating %d %d %d \n",dimp,(*p)->n[0],(*p)->n[1]);




	printf("here1\n");






	 
	    printf("here2\n");

	    //cudaMemcpy(*d_w, *w, NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);
	    cudaMemcpy(*d_wmod, *wmod, 2*(1+(((*p)->rkon)==1))*NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);
	    cudaMemcpy(*d_wd, *wd, NDERV*dimp*sizeof(real), cudaMemcpyHostToDevice);






	printf("here3\n");






	   // cudaMemcpy(*d_wnew, *wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyHostToDevice);
	    printf("here\n");
	    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
	    cudaMemcpy(*d_state, *state, sizeof(struct state), cudaMemcpyHostToDevice);
	    
	    dim3 dimBlock(16, 1);
	    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
	    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
	   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
	   

	    printf("calling initialiser\n");
	     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
	    // init_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    // init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_b);
	     init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2);


}//end of if(p->mode !=3)

     //cudaThreadSynchronize();
     

//checkErrors_i("memory allocation");

     //copy data back to cpu so we can compute and update the grid (on the cpu)

 
    //cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);

/*if((*p)->mode==3)
{
  
  int ii[3];
  ii[0]=0;
  ii[1]=0;
  ii[2]=0;
 init_user_i(*wmod,*wd,*p,ii);
}*/







    //cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);
    //setup the grid and dx values here


    //cudaMemcpy(*d_w, *w, NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);


 //initialise grid on the gpu
 //we currently don't do this to avoid use of additional memory on GPU
 /*for(dir=0; dir<NDIM; dir++)
 {
     zerotempv_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2,dir);
     cudaThreadSynchronize();     
     gridsetup_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2,dir);
     cudaThreadSynchronize();
     setupdx_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2,dir);
     cudaThreadSynchronize();
  }*/

	    printf("called initialiser\n");
	//cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);
if((*p)->mode !=3)
{
	cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);
        cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
}


//checkErrors_i("memory allocation");checkErrors_i("memory allocation");



	//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
	//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

        // printf("mod times step %f %f\n",(*p)->dt, ((*wnew)[10+16*((*p)->n[0])+((*p)->n[0])*((*p)->n[1])*b1]));



  return 0;



}




int cuupdatemod(struct params **p, struct bparams **bp,real **w, real **wnew, real **wd, struct state **state, struct params **d_p, struct bparams **d_bp,real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{
  int deviceCount;
  int dir;
 
  printf("in cuinit\n");
 // real *adb;
  real *adw, *adwnew;
  struct params *adp;
  struct bparams *adbp;
  struct state *ads;

 
 
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif  

   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
   
    printf("calling updatemod\n");
     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
    // init_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
    // init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_b);
     updatemod_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2);
     //cudaThreadSynchronize();
     


 


  return 0;
}




/*! Cartesian or polar grid. Determine x at the boundaries.
! Determine often needed combinations of x, such as dx or dvolume.
! Determine variables for axial symmetry
!
! ixe          - edge coordinate of the grid touching the boundary region
! ixf          - coordinate inside of ixe
! qx           - x with an extended index range for calculation of dx   */

int initgrid(struct params **p,   struct state **state, real **wd, struct params **d_p,  real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{
    real *ttemp2;
    int ii[NDIM];
    int ii1[3],ii2[3],ix;
    int ip,jp,kp,kpo;
    int dir,dir1,dir2;
    int ixmin,ixmax,ixe,ixf;
    real *wda=*wd;
    //real *wa=*wmod;
 int dimp=(((*p)->n[0]))*(((*p)->n[1]));

/*if(((*p)->ipe)==2)
      {
checkErrors_i("initgrid memory allocation");
}*/
 #ifdef USE_SAC_3D
 
   dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif      
    kp=0;
    //printf("called initgrid %d\n",(*p)->ipe);
    

    for(int i=0;i<3;i++)
    {
       ii1[i]=0;
       ii2[i]=0;
    }
    #ifdef USE_SAC
    ttemp2=(real *) malloc( (NTEMP2+2)*(((*p)->n[0])+2)* (((*p)->n[1])+2)*sizeof(real));
    #endif
    #ifdef USE_SAC_3D
    ttemp2=(real *)malloc((NTEMP2+2)*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2)*sizeof(real));
    #endif
    
     //cudaMemcpy(*wmod, *d_wmod, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);
     cudaMemcpy(*wd, *d_wd, NDERV*dimp*sizeof(real), cudaMemcpyDeviceToHost);
     for(dir=0;dir<NDIM;dir++)
     for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
     for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
     		     #ifdef USE_SAC_3D
                   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
                 #endif
                 {
                        ip=ii[0];
                        jp=ii[1];
         		     #ifdef USE_SAC_3D
                       kp=ii[2];
                     #endif                   
                       
	    switch(dir)
	    {
		     case 0:
	                ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui)]= 0;
		     break;
	
		     case 1:
			 ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui1)]= 0;
		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
			 ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui2)]= 0;
		     break;			
		     #endif
	     }
      }	
 

     kp=1;
     kpo=0;
     for(dir=0;dir<NDIM;dir++)
        for(ii[0]=1; ii[0]<((*p)->n[0])+1; ii[0]++)
           for(ii[1]=1; ii[1]<((*p)->n[1])+1; ii[1]++)
		#ifdef USE_SAC_3D
		   for(ii[2]=1; ii[2]<((*p)->n[2])+1; ii[2]++)
		#endif
                {
                        ip=ii[0];
                        jp=ii[1];
         		     #ifdef USE_SAC_3D
                       kp=ii[2];
                       kpo=kp;
                     #endif                   
                       
	    switch(dir)
	    {
		     case 0:
	                ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui)]= (wda[encode3_i(*p,ip-1,jp-1,kp-1,pos1)]);
		     break;
	
		     case 1:
			 ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui1)]= (wda[(encode3_i(*p,ip-1,jp-1,kp-1,pos2))]);
		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
			 ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui2)]= (wda[(encode3_i(*p,ip-1,jp-1,kp-1,pos3))]);
		     break;			
		     #endif
	     }


      }	


   /* if((*p)->ipe==3   ) 
        for(ii[1]=1; ii[1]<((*p)->n[1])+1; ii[1]++)                                                            
        for(ii[0]=1; ii[0]<((*p)->n[0])+1; ii[0]++)    
           {
                      ip=ii[0];
                        jp=ii[1];
                         printf("ii0, ii1 %d %d %16.20f %16.20f\n",ip,jp, ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui)],ttemp2[(encode3p2_i(*p,ip,jp,kpo,tmpnui1))]);

            }*/


  	
   	//update grid edges
     kp=0;
     for(dir=0;dir<NDIM;dir++)
     {
                
                       
	    switch(dir)
	    {
		     case 0:
 			    ixmax=((*p)->n[0])+1;//ixGmax1+1;
                       ixmin=((*p)->n[0])-1;//ixmin1=ixGmax1+1                      

                      #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||   ((*p)->mpiupperb[dir])==1) ixmin=((*p)->n[0])+1;//ixGmax1+1;
                      #endif

                       ixe=ixmin-1; 
                       ixf=ixe-1;


                       //upper layers
			     for(dir1=0;dir1<NDIM;dir1++)
			     {
				     for(ii[0]=ixmin; ii[0]<=ixmax; ii[0]++)
				     //for(ii[0]=ixmin; ii[0]<ixmax; ii[0]++)

				     for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
				     		 #ifdef USE_SAC_3D
						   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
						 #endif
						 {
				                        ix=ii[0];
                                                        ip=ii[0];
							jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[0]=ixe;
                                                       ii2[0]=ixf; 

 

                                                       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
						      //ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (wda[fencode3_i(*p,ii1,pos1+dir1)]);
						  }

				}



                      //lower layers

                       ixmin=0;//ixmin1=ixGmin1-1;
                      // ixmax=2;//ixmax1=ixGmin1-1 
                        ixmax=1;//ixmax1=ixGmin1-1 

                     #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||  ((*p)->mpilowerb[dir])==1) ixmax=0;
                      #endif

                  
                       ixe=ixmax+1; 
                       ixf=ixe+1;

			     for(dir1=0;dir1<NDIM;dir1++)
			     {
                                //for(ii[0]=ixmin; ii[0]<ixmax; ii[0]++)
				     for(ii[0]=ixmin; ii[0]<=ixmax; ii[0]++)
				     for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
				     		 #ifdef USE_SAC_3D
						   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
						 #endif
						 {
							ix=ip=ii[0];
							jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[0]=ixe;
                                                       ii2[0]=ixf;

 
    ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
// ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]= (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])+ (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
   // qx(ix,ixmin2:ixmax2,jdim)=(1+abs(ixe-ix))*qx(ixe,ixmin2:ixmax2,jdim)- abs(ixe-ix) *qx(ixf,ixmin2:ixmax2,jdim)


//		if((*p)->ipe==0   && ii[1]==0) 
//                         printf("ii0, ii1 %d %d %16.20f %16.20f %d %d %d %d %d  %d %d %d %d\n",ip,jp, ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui)],ttemp2[(encode3p2_i(*p,ip,jp,kpo,tmpnui1))],dir1,ixe,ixf,ixmin,ixmax,ii1[0],ii1[1],ii2[0],ii2[1]);





						  }

				}
		     break;
	
		     case 1:
                       ixmax=((*p)->n[1])+1;//ixGmax1+1; 
                       ixmin=((*p)->n[1])-1;//ixmin1=ixGmax1+1                      

                      #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||  ((*p)->mpiupperb[dir])==1) ixmin=((*p)->n[1])+1;//ixGmax1+1;
                      #endif
                     
                       ixe=ixmin-1; 
                       ixf=ixe-1;


                       //upper layers
			     for(dir1=0;dir1<NDIM;dir1++)
			     {
                 for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
				     for(ii[1]=ixmin; ii[1]<=ixmax; ii[1]++)
				     
				     		 #ifdef USE_SAC_3D
						   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
						 #endif
						 {
							ip=ii[0];
							ix=jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[1]=ixe;
                                                       ii2[1]=ixf; 





						       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
						      //ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (wda[fencode3_i(*p,ii1,pos1+dir1)]);
						  }

				}



                      //lower layers

                       ixmin=0;//ixmin1=ixGmin1-1;
                       ixmax=2;//ixmax1=ixGmin1-1 

                     #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||  ((*p)->mpilowerb[dir])==1) ixmax=0;
                      #endif
                
                       ixe=ixmax+1; 
                       ixf=ixe+1;

			     for(dir1=0;dir1<NDIM;dir1++)
			     {
			         for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)	
				     for(ii[1]=ixmin; ii[1]<=ixmax; ii[1]++)
				     		 #ifdef USE_SAC_3D
						   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
						 #endif
						 {
							ip=ii[0];
							ix=jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[1]=ixe;
                                                       ii2[1]=ixf; 



                                                        

						       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
//write(*,*) jdim,ixe,ixf,ix,ixmin1,ixmax1,ixmin2,ixmax2, qx(ixmin1:ixmax1,ix,jdim),qx(ixmin1:ixmax1,&
//                   ixe,jdim),qx(ixmin1:ixmax1,ixf,jdim)
//if((*p)->ipe==0   && ii[0]==0)
//                                                             printf("ixe, ix %d %d %d %d %d %d %16.20f %16.20f %16.20f\n",dir1,ixe,ixf,ix,ixmin,ixmax, ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)],ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))],ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
//if((*p)->ipe==0   && ii[0]==0)
//                                                             printf("ixe, ix %d %d %d %d %d %d %16.20f %16.20f %16.20f\n",dir1,ixe,ixf,ix,ixmin,ixmax, ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)],ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))],ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);


						  }

				}




		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
                       ixmax=((*p)->n[2])+1;//ixGmax1+1; 
                       ixmin=((*p)->n[2])-1;//ixmin1=ixGmax1+1                      

                      #ifdef USE_MULTIGPU
			if(((*p)->mpiupperb[dir])==1) ixmin=((*p)->n[2])+1;//ixGmax1+1;
                      #endif
                  
                       ixe=ixmin-1; 
                       ixf=ixe-1;


                       //upper layers
			     for(dir1=0;dir1<NDIM;dir1++)
			     {
                 for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
                 for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
				     
				     		 #ifdef USE_SAC_3D
						  
			        for(ii[2]=ixmin; ii[2]<=ixmax; ii[2]++)
						 #endif
						 {
							ip=ii[0];
							jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       ix=kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[2]=ixe;
                                                       ii2[2]=ixf; 
						       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
						      //ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (wda[fencode3_i(*p,ii1,pos1+dir1)]);
						  }

				}



                      //lower layers
                      //lower layers

                       ixmin=0;//ixmin1=ixGmin1-1;
                       ixmax=2;//ixmax1=ixGmin1-1 

                     #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||  ((*p)->mpilowerb[dir])==1) ixmax=0;
                      #endif
                   
                       ixe=ixmax+1; 
                       ixf=ixe+1;

			     for(dir1=0;dir1<NDIM;dir1++)
			     {
			         for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
                     for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)	
				     
				     		 #ifdef USE_SAC_3D
						   
						    for(ii[2]=ixmin; ii[2]<=ixmax; ii[2]++)
						 #endif
						 {
							ip=ii[0];
							jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       ix=kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[2]=ixe;
                                                       ii2[2]=ixf; 
						       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(fabs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
						  }

				}



		     break;			
		     #endif
	     }
      }	


kp=0;

     for(dir=0;dir<NDIM;dir++)
        for(ii[0]=0; ii[0]<((*p)->n[0]); ii[0]++)
           for(ii[1]=0; ii[1]<((*p)->n[1]); ii[1]++)
		#ifdef USE_SAC_3D
		   for(ii[2]=0; ii[2]<((*p)->n[2]); ii[2]++)
		#endif
                {
                        ip=ii[0]+1;
                        jp=ii[1]+1;
         		     #ifdef USE_SAC_3D
                       kp=ii[2]+1;
                     #endif                   
                       
	    switch(dir)
	    {
		     case 0:
	                 (wda[fencode3_i(*p,ii,pos1)])=ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui)];
                      //  if(ip==1)
                      //  printf("delx 0 %d %d %16.20f  %16.20f \n",ii[0],ii[1],wda[(encode3_i(*p,ip-1,jp-1,kp,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp,delx2))]);
		     break;
	
		     case 1:
			  (wda[(fencode3_i(*p,ii,pos2))])=ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui1)];
                 //if(ip==1)
                 //       printf("delx 1 %d %d %16.20f  %16.20f \n",ii[0],ii[1],wda[(encode3_i(*p,ip-1,jp-1,kp,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp,delx2))]);

		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
			  (wda[(fencode3_i(*p,ii,pos3))])=ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui2)];
		     break;			
		     #endif
	     }
      }	



   	//calculate dx
  kp=0;
  kpo=0;

                   for(dir=0;dir<NDIM;dir++)
                 {

  for(ii[0]=1; ii[0]<((*p)->n[0])+1; ii[0]++)
     for(ii[1]=1; ii[1]<((*p)->n[1])+1; ii[1]++)
  //for(ii[0]=0; ii[0]<((*p)->n[0]); ii[0]++)
  //   for(ii[1]=0; ii[1]<((*p)->n[1]); ii[1]++)

     		     #ifdef USE_SAC_3D
                   for(ii[2]=1; ii[2]<((*p)->n[2])+1; ii[2]++)
                 #endif
{

                        ip=ii[0];
                        jp=ii[1];
         		     #ifdef USE_SAC_3D
                       
                       kp=ii[2];
                        kpo=kp-1;
                     #endif                   
                       
	    switch(dir)
	    {
		     case 0:
	               (wda[(encode3_i(*p,ip-1,jp-1,kpo,delx1))])=/*(*p)->dx[0];//*/0.5*(ttemp2[encode3p2_i(*p,ip+1,jp,kp,tmpnui)]-ttemp2[encode3p2_i(*p,ip-1,jp,kp,tmpnui)]);
                      //(wda[(encode3_i(*p,ip-1,jp-1,kpo,delx1))])=/*(*p)->dx[0];//*/0.5*(ttemp2[encode3p2_i(*p,ip+1,jp,kp,tmpnui)]-ttemp2[encode3p2_i(*p,ip-1,jp,kp,tmpnui)]);
	              //  if(ip==128  && jp==128 && kp==128)
                      //  printf("delx 0 %d %d %d %16.20f  %16.20f   %16.20f \n",ii[0]-1,ii[1]-1,ii[2]-1,wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx2))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx3))]);
		     break;
	
		     case 1:
			(wda[(encode3_i(*p,ip-1,jp-1,kpo,delx2))])=/*(*p)->dx[1];//*/0.5*(ttemp2[encode3p2_i(*p,ip,jp+1,kp,tmpnui1)]-ttemp2[encode3p2_i(*p,ip,jp-1,kp,tmpnui1)]);
			//(wda[(encode3_i(*p,ip-1,jp-1,kpo,delx2))])=/*(*p)->dx[1];//*/0.5*(ttemp2[encode3p2_i(*p,ip,jp+1,kp,tmpnui1)]-ttemp2[encode3p2_i(*p,ip,jp-1,kp,tmpnui1)]);

	               // if(ip==128  && jp==128 && kp==128)
                       //   printf("delx 1 %d %d %d %16.20f  %16.20f   %16.20f \n",ii[0]-1,ii[1]-1,ii[2]-1,wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx2))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx3))]);

		        //printf("delx2 %d %d %g ",ii[0],ii[1],wda[(fencode3_i(*p,ii,delx2))]);
		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
			(wda[(encode3_i(*p,ip-1,jp-1,kpo,delx3))])=0.5*(ttemp2[encode3p2_i(*p,ip,jp,kp+1,tmpnui2)]-ttemp2[encode3p2_i(*p,ip,jp,kp-1,tmpnui2)]);
	              //  if(ip==128  && jp==128 && kp==128)
                      //  printf("delx 2 %d %d %d %16.20f  %16.20f   %16.20f \n",ii[0]-1,ii[1]-1,ii[2]-1,wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx2))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx3))]);

		     break;			
		     #endif
	     }
      }
  printf("\n");
}


printf("dx=%g dy=%g\n",(*p)->dx[0], (*p)->dx[1] );




     kp=0;

   // if((*p)->ipe==3)
    // for(dir=0;dir<NDIM;dir++)
//for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
        //for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
           
            // {

              //          ip=ii[0];
              //          jp=ii[1];
                //if(ii[0]==0)
              //  printf("delx %d %d %16.20f  %16.20f  %16.20f  %16.20f \n",ii[0],ii[1],wda[(fencode3_i(*p,ii,pos1))],wda[(fencode3_i(*p,ii,pos2))],wda[(fencode3_i(*p,ii,delx1))],wda[(fencode3_i(*p,ii,delx2))]);
//printf("ttemp2 %d %d %16.20f  %16.20f  \n",ii[0],ii[1],ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui)],ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui1)]);
            //  }

  //  cudaMemcpy(*d_w, *w, NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);
  //  cudaMemcpy(*d_wd, *wd, NDERV*dimp*sizeof(real), cudaMemcpyHostToDevice);

 free(ttemp2);

 

    //cudaMemcpy(*d_wmod, *wmod, NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_wd, *wd, NDERV*dimp*sizeof(real), cudaMemcpyHostToDevice);


    //  
     



  

   
  return 0;



}


#ifdef USE_MULTIGPU

//prepare data buffers used to copy data between gpu and cpu
//this will update only the ghost cells transferred between the CPU's


int cuinitmgpubuffers(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,   real **gmpiw0, real **gmpiwmod0,   real **gmpiw1, real **gmpiwmod1,   real **gmpiw2, real **gmpiwmod2, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2)
{

  int szw,  szvisc0,szvisc1,szvisc2,szw0,szw1,szw2;
  #ifdef USE_SAC
  //real *dt;
  
  szw=4*(  ((*p)->n[1])  +  ((*p)->n[0])   );
  szw0=4*NDERV*(  ((*p)->n[1])     );
  szw1=4*NDERV*(  ((*p)->n[0])     );

  szvisc0=4*(  (((*p)->n[1])+2 )   );
  szvisc1=4*(    (((*p)->n[0]) +2 )  );

 //dt=(real *)calloc( NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2),sizeof(real));

  #endif
  #ifdef USE_SAC_3D
  
  szw=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])  +  ((*p)->n[0])*((*p)->n[2])  +  ((*p)->n[0])*((*p)->n[1])  );
  szw0=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NDERV*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NDERV*(    ((*p)->n[0])*((*p)->n[1])  );



  szvisc0=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)  ); 
  szvisc1=4*(   (((*p)->n[0])+2)*(((*p)->n[2])+2)    );    
  szvisc2=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)   );    

   
  //dt=(real *)calloc( NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2),sizeof(real));
  #endif






  //temp2=&dt;
  //gmpiwmod=(real **)malloc(szw*sizeof(real));
  //gmpiw=(real **)malloc(szw*sizeof(real));

  //gmpiwmod0=(real **)malloc(szw0*sizeof(real));
  //gmpiw0=(real **)malloc(szw0*sizeof(real));
  //gmpiwmod1=(real **)malloc(szw1*sizeof(real));
  //gmpiw1=(real **)malloc(szw1*sizeof(real));

  #ifdef USE_SAC_3D
	//  gmpiwmod2=(real **)malloc(szw2*sizeof(real));
	//  gmpiw2=(real **)malloc(szw2*sizeof(real));
  #endif

  //gmpivisc=(real **)malloc(szvisc*sizeof(real));
	//  cudaMalloc((void**)d_gmpiwmod, szw*sizeof(real));
	//  cudaMalloc((void**)d_gmpiw, szw*sizeof(real));


  	  cudaMalloc((void**)d_gmpiwmod0, szw0*sizeof(real));
	  cudaMalloc((void**)d_gmpiw0, szw0*sizeof(real));
	  cudaMalloc((void**)d_gmpiwmod1, szw1*sizeof(real));
	  cudaMalloc((void**)d_gmpiw1, szw1*sizeof(real));

  #ifdef USE_SAC_3D  
	  cudaMalloc((void**)d_gmpiwmod2, szw2*sizeof(real));
	  cudaMalloc((void**)d_gmpiw2, szw2*sizeof(real));
          cudaMalloc((void**)d_gmpivisc2, szvisc2*sizeof(real));
  #else

          cudaMalloc((void**)d_gmpivisc2, sizeof(real));
  #endif
          cudaMalloc((void**)d_gmpivisc0, szvisc0*sizeof(real));
          cudaMalloc((void**)d_gmpivisc1, szvisc1*sizeof(real));
  return 0;
}

//copy gpu memory data to mpi send buffer for w and wmod
//just update the edges of w and wmod with values copied from gmpiw, gmpiwmod and gmpivisc
int cucopywtompiwmod(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir)
{
     int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;
     int szw0,szw1,szw2;

     int szbuf;
     int dimp=(((*p)->n[0]))*(((*p)->n[1]));
     
     
     i3=0;
     #ifdef USE_SAC_3D  
       dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
     #endif 
     int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif


  #ifdef USE_SAC
  
  szw0=4*NVAR*(  ((*p)->n[1])     );
  szw1=4*NVAR*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NVAR*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NVAR*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NVAR*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif

    //real **d_tgmpiw0;
    //real **tgmpiw0=(real **)malloc(szw0*sizeof(real));
    //cudaMalloc((void**)d_tgmpiw0, szw0*sizeof(real));
    // for(var=0; var<NVAR; var++)
    //   for(dim=0;dim<NDIM;dim++)
     gputompiwmod_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod,*d_gmpiw0,*d_gmpiwmod0,*d_gmpiw1,*d_gmpiwmod1,*d_gmpiw2,*d_gmpiwmod2,order,idir);

     
     cudaThreadSynchronize();

if(idir==0)
{
     cudaMemcpy(*gmpiwmod0, *d_gmpiwmod0, szw0*sizeof(real), cudaMemcpyDeviceToHost);

// cudaThreadSynchronize();

   //   cudaMemcpy(*gmpiw0, *d_gmpiw0, szw0*sizeof(real), cudaMemcpyDeviceToHost);




}

if(idir==1)
{
     cudaMemcpy(*gmpiwmod1, *d_gmpiwmod1, szw1*sizeof(real), cudaMemcpyDeviceToHost);
   //  cudaMemcpy(*gmpiw1, *d_gmpiw1, szw1*sizeof(real), cudaMemcpyDeviceToHost);
}
 

//struct params *tp;     







    

   #ifdef USE_SAC3D
if(idir==2)
{
     cudaMemcpy(*gmpiwmod2, *d_gmpiwmod2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*gmpiw2, *d_gmpiw2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
}
   #endif 

cudaThreadSynchronize();
}


//copy gpu memory data to mpi send buffer for w and wmod
//just update the edges of w and wmod with values copied from gmpiw, gmpiwmod and gmpivisc
int cucopywtompiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir)
{
     int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;
     int szw0,szw1,szw2;

     int szbuf;
     int dimp=(((*p)->n[0]))*(((*p)->n[1]));
     
     
     i3=0;
     #ifdef USE_SAC_3D  
       dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
     #endif 
     int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif


  #ifdef USE_SAC
  
  szw0=4*NVAR*(  ((*p)->n[1])     );
  szw1=4*NVAR*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NVAR*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NVAR*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NVAR*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif

    //real **d_tgmpiw0;
    //real **tgmpiw0=(real **)malloc(szw0*sizeof(real));
    //cudaMalloc((void**)d_tgmpiw0, szw0*sizeof(real));
    // for(var=0; var<NVAR; var++)
    //   for(dim=0;dim<NDIM;dim++)
     gputompiw_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod,*d_gmpiw0,*d_gmpiwmod0,*d_gmpiw1,*d_gmpiwmod1,*d_gmpiw2,*d_gmpiwmod2,order,idir);

     
     cudaThreadSynchronize();

if(idir==0)
{
     //cudaMemcpy(*gmpiwmod0, *d_gmpiwmod0, szw0*sizeof(real), cudaMemcpyDeviceToHost);


      cudaMemcpy(*gmpiw0, *d_gmpiw0, szw0*sizeof(real), cudaMemcpyDeviceToHost);
}

if(idir==1)
{
     //cudaMemcpy(*gmpiwmod1, *d_gmpiwmod1, szw1*sizeof(real), cudaMemcpyDeviceToHost);
     cudaMemcpy(*gmpiw1, *d_gmpiw1, szw1*sizeof(real), cudaMemcpyDeviceToHost);
}
      
    

   #ifdef USE_SAC3D
if(idir==2)
{
    // cudaMemcpy(*gmpiwmod2, *d_gmpiwmod2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
     cudaMemcpy(*gmpiw2, *d_gmpiw2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
}
   #endif 

cudaThreadSynchronize();


 /*if(((*p)->ipe)==3  && ((*p)->it)==2)
{


       printf("%d %d \n",szw0,szw1);

        for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=0;
                       bound=0;
                       var=0;
                            printf(" %d %d %d %lg %d \n",i1,i2,bound,(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)],encodempiw1(*p,i1,i2,i3,var,bound));                                        
                     ;//  (*wmod)[fencode3_i(*p,ii,var)]=(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)];              
                     ;//  (*w)[fencode3_i(*p,ii,var)]=(*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)];



                  }
         ;// printf("\n");
}*/









//encodempiw1(p,i,j,k,var,bound)]



   //free(tgmpiw0);
   //cudaFree(*d_tgmpiw0);
//gmpiw behaving OK but cannot display or access any of the gmpiwmod variables!
//printf("%f\n",(*gmpiwmod)[0]);
     
//encodempiw (struct params *dp,int ix, int iy, int iz, int field,int bound,int dim)
     //copy data to correct area in w and wmod
   /*  for(var=0; var<NVAR; var++)
       for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<4;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC3D
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;                                                                     
                       (*wmod)[fencode3_i(*p,ii,var)]=(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)];              
                       (*w)[fencode3_i(*p,ii,var)]=(*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)];
                  }
            #else
         ii[2]=0;
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                      


		
                       (*wmod)[fencode3_i(*p,ii,var)]=(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)];  
                      //if(encodempiw(*p,i1,i2,i3,var,bound,dim)<10)  
                      if(var==5 && ((*p)->ipe)==0)                        
                        printf(" %d %d %d %d actual %d  mpi data%d %g\n",i1,i2,bound,dim,var,encodempiw0(*p,i1,i2,i3,var,bound),(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]);

                     ;// if(encodempiw(*p,i1,i2,i3,var,bound,dim)<10239 )
                       (*w)[fencode3_i(*p,ii,var)]=(*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)];
                                                                     
                      // *(wmod+encode3_i(*p,ii,var))=*(gmpiwmod0+encodempiw0(*p,i1,i2,i3,var,bound));              
                      // (*w)[encode3_i(*p,ii,var)]=(*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)];

                      /* if(var==4  && ((*p)->ipe)==1)
                       {
				(*wmod)[fencode3_i(*p,ii,var)]=0.5;
				(*w)[fencode3_i(*p,ii,var)]=0.5;
                       }*/


              /*    }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC3D
         i2=bound*(bound<2)+(((*p)->n[1])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;                                                                     
                       (*wmod)[fencode3_i(*p,ii,var)]=(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)];              
                       (*w)[fencode3_i(*p,ii,var)]=(*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)];
                  }

            #else
         ii[2]=0;
         i2=bound*(bound<2)+(   ((*p)->n[1])-(bound-1)   )*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                                                                     
                     ;//  (*wmod)[fencode3_i(*p,ii,var)]=(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)];              
                     ;//  (*w)[fencode3_i(*p,ii,var)]=(*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)];



                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC3D
                       case 2:
         i3=bound*(bound<2)+( ((*p)->n[2])-(bound-1) )*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;                                                                     
                       (*wmod)[fencode3_i(*p,ii,var)]=(*gmpiwmod2)[encodempiw2(*p,i1,i2,i3,var,bound)];              
                       (*w)[fencode3_i(*p,ii,var)]=(*gmpiw2)[encodempiw2(*p,i1,i2,i3,var,bound)];
                  }                            
                       break;                       
            #endif             
             }
                                     
         }    */

}




int cucopywdtompiwd(struct params **p,real **wd,    real **gmpiw0,    real **gmpiw1,    real **gmpiw2, struct params **d_p  ,real **d_wd,   real **d_gmpiw0,   real **d_gmpiw1,   real **d_gmpiw2, int order, int idir)
{
     int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;
     int szw0,szw1,szw2;

     int szbuf;
     int dimp=(((*p)->n[0]))*(((*p)->n[1]));
     
     
     i3=0;
     #ifdef USE_SAC_3D  
       dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
     #endif 
     int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif


  #ifdef USE_SAC
  
  szw0=4*NDERV*(  ((*p)->n[1])     );
  szw1=4*NDERV*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NDERV*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NDERV*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif

    //real **d_tgmpiw0;
    //real **tgmpiw0=(real **)malloc(szw0*sizeof(real));
    //cudaMalloc((void**)d_tgmpiw0, szw0*sizeof(real));
    // for(var=0; var<NVAR; var++)
    //   for(dim=0;dim<NDIM;dim++)
     gputompiwd_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_wd,*d_gmpiw0,*d_gmpiw1,*d_gmpiw2,order, idir);

     
     cudaThreadSynchronize();
if(idir==0)
      cudaMemcpy(*gmpiw0, *d_gmpiw0, szw0*sizeof(real), cudaMemcpyDeviceToHost);

if(idir==1)
     cudaMemcpy(*gmpiw1, *d_gmpiw1, szw1*sizeof(real), cudaMemcpyDeviceToHost);

      
    

   #ifdef USE_SAC3D
if(idir==2)
      cudaMemcpy(*gmpiw2, *d_gmpiw2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
   #endif 

cudaThreadSynchronize();


 /*if(((*p)->ipe)==3  && ((*p)->it)==2)
{


       printf("%d %d \n",szw0,szw1);

        for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=0;
                       bound=0;
                       var=0;
                            printf(" %d %d %d %lg %d \n",i1,i2,bound,(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)],encodempiw1(*p,i1,i2,i3,var,bound));                                        
                     ;//  (*wmod)[fencode3_i(*p,ii,var)]=(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)];              
                     ;//  (*w)[fencode3_i(*p,ii,var)]=(*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)];



                  }
         ;// printf("\n");
}*/





}




//copy mpi recv buffer to gpu memory     
int cucopywfrommpiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir)
{
       int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;     
       int szbuf;
     int szw0,szw1,szw2;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D  
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif      
     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif
        int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


  #ifdef USE_SAC
  
  szw0=4*NVAR*(  ((*p)->n[1])     );
  szw1=4*NVAR*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NVAR*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NVAR*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NVAR*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif




      //copy data from w and wmod to correct gmpiw and gmpiwmod

//encodempiw (struct params *dp,int ix, int iy, int iz, int field,int bound,int dim)
     //copy data to correct area in w and wmod
   /*  for(var=0; var<NVAR; var++)
       for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<4;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC3D
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;                                                                     
                       (*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];              
                       (*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];
                  }
            #else
         ii[2]=0;
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       (*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];  

                 //   if(var==4 && ((*p)->ipe)==1) 
                                 
                //        printf(" %d %d %d %d actual %d  mpi data%d %g %g\n",i1,i2,bound,dim,var,encodempiw0(*p,i1,i2,i3,var,bound),(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)],  (*wmod)[fencode3_i(*p,ii,var)] );
            
                       (*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];

                      //  if(var==4  && ((*p)->ipe)==1)
                     //  {
			//	(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=0.5;
			//	(*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=0.5;
                      // }

                  }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC3D
         i2=bound*(bound<2)+(((*p)->n[1])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;  

                       (*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];              
                       (*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];

                  }

            #else
         ii[2]=0;
         i2=bound*(bound<2)+(   ((*p)->n[1])-(bound-1)   )*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                      (*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];              
                       (*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];      

                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC3D
                       case 2:
         i3=bound*(bound<2)+(((*p)->n[2])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3; 

                      (*gmpiwmod2)[encodempiw2(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];              
                       (*gmpiw2)[encodempiw2(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];      
                    }                            
                       break;                       
            #endif             
             }
                                     
         }    */

//encodempiw (struct params *dp,int ix, int iy, int iz, int field,int bound,int dim)
     //copy data to correct area in w and wmod
     /*for(var=0; var<NVAR; var++)
       for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<4;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC3D
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;     

                      (*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];              
                       (*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];      
  
                  }
            #else
         ii[2]=0;
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;

                      (*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];              
                       (*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];      
                  }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC3D
         i2=bound*(bound<2)+(((*p)->n[1])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3; 

                      (*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];              
                       (*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];      
                   }

            #else
         ii[2]=0;
         i2=bound*(bound<2)+(   ((*p)->n[1])-(bound-1)   )*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;


                      (*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];              
                       (*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];      
                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC3D
                       case 2:
         i3=bound*(bound<2)+(((*p)->n[2])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3; 


                      (*gmpiwmod2)[encodempiw2(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_i(*p,ii,var)];              
                       (*gmpiw2)[encodempiw2(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_i(*p,ii,var)];      
                   }                            
                       break;                       
            #endif             
             }
                                     
         }  */ 



if(idir==0)
{
   	 cudaMemcpy(*d_gmpiw0, *gmpiw0, szw0*sizeof(real), cudaMemcpyHostToDevice);     
   	 //cudaMemcpy(*d_gmpiwmod0, *gmpiwmod0, szw0*sizeof(real), cudaMemcpyHostToDevice); 
}

if(idir==1)
{
   	 cudaMemcpy(*d_gmpiw1, *gmpiw1, szw1*sizeof(real), cudaMemcpyHostToDevice);     
   	 //cudaMemcpy(*d_gmpiwmod1, *gmpiwmod1, szw1*sizeof(real), cudaMemcpyHostToDevice);     
}
    
            #ifdef USE_SAC3D
if(idir==2)
{
   	      cudaMemcpy(*d_gmpiw2, *gmpiw2, szw2*sizeof(real), cudaMemcpyHostToDevice);     
   	     // cudaMemcpy(*d_gmpiwmod2, *gmpiwmod0, szw2*sizeof(real), cudaMemcpyHostToDevice); 
}    
         #endif

    //printf("call mpiwtogpu\n");

     mpiwtogpu_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod,*d_gmpiw0,*d_gmpiwmod0,*d_gmpiw1,*d_gmpiwmod1,*d_gmpiw2,*d_gmpiwmod2,idir);
     cudaThreadSynchronize();
}

int cucopywmodfrommpiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir)
{
       int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;     
       int szbuf;
     int szw0,szw1,szw2;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));



real *tgmpiwmod0=*gmpiwmod0;
real *tgmpiwmod1=*gmpiwmod1;
 /*  if(idir==0 && (*p)->ipe==0     &&  idir==0 )
    {
        printf("ipe2 mpiw0 after bound \n");
        for(int j=0; j<4;j++)
         for(int i=0; i<((*p)->n[1]);i++)           
             printf("%d %d %lg %lg\n",i,j, (tgmpiwmod0[4*rho*((*p)->n[0]) +i+j*((*p)->n[0])]), (tgmpiwmod1[4*rho*((*p)->n[0]) +i+j*((*p)->n[0])]));
         printf("\n");
     }*/


   
 #ifdef USE_SAC_3D  
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif      
     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif
        int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


  #ifdef USE_SAC
  
  szw0=4*NVAR*(  ((*p)->n[1])     );
  szw1=4*NVAR*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NVAR*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NVAR*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NVAR*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif




      //copy data from w and wmod to correct gmpiw and gmpiwmod




if(idir==0)
{



		    i3=0;
                   //   int bound,iside,n=0;
                     //for(iside=0;iside<2;iside++)
		    //for(int i1=0;i1<=1;i1++)
		      //for(int i2=0;i2<(*p)->n[1];i2++)
		      //{
			//iside=0;
                        // bound=i1+2*(iside>0);
			 
			//if(((*p)->ipe==0) /*&&  (*p)->it != -1     && iside==1 && (100*(p->ipe)+10*dim+iside)==101*/ )
			//{
                          //  printf("tini %d %d %d %lg  \n",bound,i2,i1,tgmpiwmod0[encodempiw0 (*p,i1, i2, i3, 0,bound)]);
                            // printf(" %d %d %d   \n",bound,i2,i1);

			//}
                   //   n++;
                    //}





   	 cudaMemcpy(*d_gmpiwmod0, *gmpiwmod0, szw0*sizeof(real), cudaMemcpyHostToDevice); 
}

if(idir==1)
{
   	 cudaMemcpy(*d_gmpiwmod1, *gmpiwmod1, szw1*sizeof(real), cudaMemcpyHostToDevice);     
}
    
            #ifdef USE_SAC3D
if(idir==2)
{
   	      cudaMemcpy(*d_gmpiwmod2, *gmpiwmod0, szw2*sizeof(real), cudaMemcpyHostToDevice); 
}    
         #endif

    //printf("call mpiwtogpu\n");

     mpiwmodtogpu_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod,*d_gmpiw0,*d_gmpiwmod0,*d_gmpiw1,*d_gmpiwmod1,*d_gmpiw2,*d_gmpiwmod2,idir,order);
     cudaThreadSynchronize();
}


//copy mpi recv buffer to gpu memory     
int cucopywdfrommpiwd(struct params **p,real **wd,     real **gmpiw0,     real **gmpiw1,     real **gmpiw2,  struct params **d_p  ,real **d_wd,    real **d_gmpiw0,   real **d_gmpiw1,   real **d_gmpiw2,  int order, int idir)
{
       int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;     
       int szbuf;
     int szw0,szw1,szw2;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D  
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif      
     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif
        int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


  #ifdef USE_SAC
  
  szw0=4*NDERV*(  ((*p)->n[1])     );
  szw1=4*NDERV*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NDERV*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NDERV*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif




       if(idir==0)
   	 cudaMemcpy(*d_gmpiw0, *gmpiw0, szw0*sizeof(real), cudaMemcpyHostToDevice);     
 
	if(idir==1)
   	 cudaMemcpy(*d_gmpiw1, *gmpiw1, szw1*sizeof(real), cudaMemcpyHostToDevice);     

    
            #ifdef USE_SAC3D
     if(idir==2)
   	      cudaMemcpy(*d_gmpiw2, *gmpiw2, szw2*sizeof(real), cudaMemcpyHostToDevice);     
         #endif

    //printf("call mpiwtogpu\n");

     mpiwdtogpu_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,0,0,*d_wd,*d_gmpiw0,*d_gmpiw1,*d_gmpiw2,idir);
     cudaThreadSynchronize();
}



//copy gpu memory data to mpi send buffer for w and wmod
//just update the edges of w and wmod with values copied from gmpiw, gmpiwmod and gmpivisc
int cucopytompivisc(struct params **p,real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2)
{


     int szbuf,szbuf0,szbuf1,szbuf2;
     int dim,bound,var=0;
     int i1,i2,i3;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
             int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
	  szbuf0=4*(  (((*p)->n[1])+2 )   );
	  szbuf1=4*(    (((*p)->n[0]) +2 )  );


     #ifdef USE_SAC3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
  szbuf0=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)  ); 
  szbuf1=4*(   (((*p)->n[0])+2)*(((*p)->n[2])+2)    );    
  szbuf2=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)   );    


     #endif
     gputompivisc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_wtemp2,*d_gmpivisc0,*d_gmpivisc1,*d_gmpivisc2);
     cudaThreadSynchronize();
     cudaMemcpy(*gmpivisc0, *d_gmpivisc0, szbuf0*sizeof(real), cudaMemcpyDeviceToHost);
     cudaMemcpy(*gmpivisc1, *d_gmpivisc1, szbuf1*sizeof(real), cudaMemcpyDeviceToHost);
     #ifdef USE_SAC3D
     	cudaMemcpy(*gmpivisc2, *d_gmpivisc2, szbuf2*sizeof(real), cudaMemcpyDeviceToHost);
     #endif
     //copy data to correct area in temp2
//encodempiw (struct params *dp,int ix, int iy, int iz, int field,int bound,int dim)
     //copy data to correct area in w and wmod
     /*  for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<2;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC3D
         i1=bound*(((*p)->n[0])+1);
         for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  for(i3=1;i3<(((*p)->n[2])+2);i3++ )
                  {     
                        
          //i1=(p->n[0])+1;
         
          //temp2[encode3p2_sacmpi (p,i1, i2, i3, tmpnui)]=gmpitgtbufferr[0][i2+i3*((p->n[1])+2)];
          //temp2[encode3p2_sacmpi (p,0, i2, i3, tmpnui)]=gmpitgtbufferl[0][i2+i3*((p->n[1])+2)];
         
                       (*temp2)[encode3p2_i(*p,i1,i2,i3,var)]=(*gmpivisc0)[encodempivisc0(*p,i1,i2,i3,bound,dim)];
                  }
            #else
         i3=0;
         i1=bound*(((*p)->n[0])+1);
                  for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  {
                       (*temp2)[encode3p2_i(*p,i1,i2,i3,var)]=(*gmpivisc0)[encodempivisc0(*p,i1,i2,i3,bound,dim)];
                  }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC3D
         i2=bound*(((*p)->n[1])+1);
         for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  for(i3=1;i3<(((*p)->n[2])+2);i3++ )
                  {
                       (*temp2)[encode3p2_i(*p,i1,i2,i3,var)]=(*gmpivisc1)[encodempivisc1(*p,i1,i2,i3,bound,dim)];
                  }

            #else
         i3=0;
         i2=bound*(((*p)->n[1])+1);
                  for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  {                                                       
                       (*temp2)[encode3p2_i(*p,i1,i2,i3,var)]=(*gmpivisc1)[encodempivisc1(*p,i1,i2,i3,bound,dim)];
                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC3D
                       case 2:
                  i3=bound*(((*p)->n[2])+1);
        for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  {                                                           
                       (*temp2)[encode3p2_i(*p,i1,i2,i3,var)]=(*gmpivisc2)[encodempivisc2(*p,i1,i2,i3,bound,dim)];
                  }                            
                       break;                       
            #endif             
             }
                                     
         }    */

}

//copy mpi recv buffer to gpu memory     
int cucopyfrommpivisc(struct params **p,real **temp2,real **gmpivisc0,real **gmpivisc1,real **gmpivisc2,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2)
{
      int dim,bound,var=0;
     int i1,i2,i3;      
 
     int szbuf,szbuf0,szbuf1,szbuf2;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D  
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

        int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

     
     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     
 	  szbuf0=4*(  (((*p)->n[1])+2 )   );
	  szbuf1=4*(    (((*p)->n[0]) +2 )  );

     #ifdef USE_SAC3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );

     
  szbuf0=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)  ); 
  szbuf1=4*(   (((*p)->n[0])+2)*(((*p)->n[2])+2)    );    
  szbuf2=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)   ); 

     #endif

      //copy data from temp2 to gmpivisc
        /*     for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<2;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC3D
         i1=bound*(((*p)->n[0])+1);
         for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  for(i3=1;i3<(((*p)->n[2])+2);i3++ )
                  {     
                        
          //i1=(p->n[0])+1;
         
          //temp2[encode3p2_sacmpi (p,i1, i2, i3, tmpnui)]=gmpitgtbufferr[0][i2+i3*((p->n[1])+2)];
          //temp2[encode3p2_sacmpi (p,0, i2, i3, tmpnui)]=gmpitgtbufferl[0][i2+i3*((p->n[1])+2)];
         
                       (*gmpivisc0)[encodempivisc0(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_i(*p,i1,i2,i3,var)];
                  }
            #else
         i3=0;
         i1=bound*(((*p)->n[0])+1);
                  for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  {
                       (*gmpivisc0)[encodempivisc0(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_i(*p,i1,i2,i3,var)];
                  }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC3D
         i2=bound*(((*p)->n[1])+1);
         for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  for(i3=1;i3<(((*p)->n[2])+2);i3++ )
                  {
                       (*gmpivisc1)[encodempivisc1(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_i(*p,i1,i2,i3,var)];
                  }

            #else
         i3=0;
         i2=bound*(((*p)->n[1])+1);
                  for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  {
                                                                     
                       (*gmpivisc1)[encodempivisc1(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_i(*p,i1,i2,i3,var)];
                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC3D
                       case 2:
                  i3=bound*(((*p)->n[2])+1);
        for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  {
                                                              
                       (*gmpivisc2)[encodempivisc2(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_i(*p,i1,i2,i3,var)];
                  }                            
                       break;                       
            #endif             
             }
                                     
         } */   


   	 cudaMemcpy(*d_gmpivisc0, *gmpivisc0, szbuf0*sizeof(real), cudaMemcpyHostToDevice);     
   	 cudaMemcpy(*d_gmpivisc1, *gmpivisc1, szbuf1*sizeof(real), cudaMemcpyHostToDevice);
       #ifdef USE_SAC3D    
   	 cudaMemcpy(*d_gmpivisc2, *gmpivisc2, szbuf2*sizeof(real), cudaMemcpyHostToDevice); 
       #endif    

     mpivisctogpu_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_wtemp2,*d_gmpivisc0,*d_gmpivisc1,*d_gmpivisc2);
     cudaThreadSynchronize();
}


#endif



