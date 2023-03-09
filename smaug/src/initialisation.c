

#include "../include/initialisation.h"






int encode3_in(Params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
  #else
    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
  #endif
}

void initconfig(Params *k, Meta *md, real *w, real *wd)
{




	int i1,j1,k1;
        int ni=k->n[0];
        int nj=k->n[1];
        unsigned long int ilv;
        printf("%d %d\n",ni,nj);
        for(i1=0; i1<(k->n[0]) ;i1++)
	  for(j1=0; j1<(k->n[1]) ;j1++)
          {

                    wd[encode3_in(k,i1,j1,k1,pos1)]=(k->xmin[0])+((real)i1)*((k->xmax[0])- (k->xmin[0])  )/ni;
                    wd[encode3_in(k,i1,j1,k1,delx1)]=((k->xmax[0])- (k->xmin[0])  )/ni;
                    wd[encode3_in(k,i1,j1,k1,pos2)]=(k->xmin[1])+((real)j1)*((k->xmax[1])- (k->xmin[1])  )/nj;
                    wd[encode3_in(k,i1,j1,k1,delx2)]=((k->xmax[1])- (k->xmin[1])  )/nj;

		    #ifdef USE_SAC3D
		            wd[encode3_in(k,i1,j1,k1,pos3)]=(k->xmin[2])+((real)k1)*((k->xmax[2])-(k->xmin[2]))/nk;
		            wd[encode3_in(k,i1,j1,k1,delx3)]=((k->xmax[2])- (k->xmin[2])  )/nk;
                    #endif

                    for(int f=rho; f<NVAR; f++)
                    //for(int f=rho; f<=b2; f++)
                    {
                    ;//ilv=j1*ni+i1+(ni*nj*f);
                     ilv=j1*ni+i1+(ni*nj*f);
                    w[ilv]=0.0;
                    switch(f)
		            {
		              case rho:
		            	w[ilv]=1.0;
			      break;
		              case mom1:
		            	w[ilv]=0.01;
			      break;
		              case mom2:
		            	w[ilv]=0.01;
			      break;
		              //case mom3:
		            	//w[j1*ni+i1+(ni*nj*f)]=0.0;
			      //break;
		              case energy:
		            	w[ilv]=0.0;
			      break;
		              case b1:
		            	w[ilv]=0.0;
			      break;
		              case b2:
		            	w[ilv]=0.0;
			      break;
		              //case b3:
		            //	w[j1*ni+i1+(ni*nj*f)]=0.0;
			     // break;
		            }; //end of switch to check for field

			}//end of loop over f

          }//end of loop over j and i
          printf("ilv %ld\n", ilv);


}

int encode3_uin(Params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
  #else
    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
  #endif
}

int fencode3_uin (Params *dp,int *ii, int field) {


#ifdef USE_SAC_3D
   return (ii[2]*((dp)->n[0])*((dp)->n[1])  + ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
#else
   return ( ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])));
#endif

}


#ifdef USE_SAC_3D
real grad3dngen_uin(real ***wmod, real *wd,struct params *p,int *ii,int dir)
#else
real grad3dngen_uin(real **wmod, real *wd,Params *p,int *ii,int dir)
#endif
{


 real grad=0;



 switch(dir)
 {
   case 0:
#ifdef USE_SAC_3D

if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[ii[0]+1][ii[1]][ii[2]]-8*wmod[ii[0]-1][ii[1]][ii[2]]+wmod[ii[0]-2][ii[1]][ii[2]]-wmod[ii[0]+2][ii[1]][ii[2]])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx1)]))    );




  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
  #endif

#else


if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[ii[0]+1][ii[1]]-8*wmod[ii[0]-1][ii[1]]+wmod[ii[0]-2][ii[1]]-wmod[ii[0]+2][ii[1]])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx1)]))    );



  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
  #endif
#endif


   break;

   case 1:


#ifdef USE_SAC_3D


if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[ii[0]][ii[1]+1][ii[2]]-8*wmod[ii[0]][ii[1]-1][ii[2]]+wmod[ii[0]][ii[1]-2][ii[2]]-wmod[ii[0]][ii[1]+2][ii[2]])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx2)]))    );


  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
  #endif
#else


if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[ii[0]][ii[1]+1]-8*wmod[ii[0]][ii[1]-1]+wmod[ii[0]][ii[1]-2]-wmod[ii[0]][ii[1]+2])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx2)]))    );



  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
  #endif
#endif


   break;


   case 2:
#ifdef USE_SAC_3D

if( ii[2] >1 &&  ii[2]<((p->n[2])-2))
	grad=(  ( ((8*wmod[ii[0]][ii[1]][ii[2]+1]-8*wmod[ii[0]][ii[1]][ii[2]-1]+wmod[ii[0]][ii[1]][ii[2]-2]-wmod[ii[0]][ii[1]][ii[2]+2])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx3)]))    );




  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[2]==(p->n[2])-3) || (ii[2]==(p->n[2])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[2]==2 || ii[2]==3  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
  #endif
#endif

   break;


}



 return grad;


}



real grad3dn_uin(real *wmod, real *wd,struct params *p,int *ii,int field,int dir)
{


 real grad=0;



 switch(dir)
 {
   case 0:

#ifdef USE_SAC_3D
  #ifdef USE_DORDER3
 if(ii[0]>2 && ii[0]<((p->n[0])-3) )
  grad=(  ( ((3*wmod[encode3_uin(p,ii[0]+1,ii[1],ii[2],field)]-3*wmod[encode3_uin(p,ii[0]-1,ii[1],ii[2],field)]+3.0*(wmod[encode3_uin(p,ii[0]-2,ii[1],ii[2],field)]-wmod[encode3_uin(p,ii[0]+2,ii[1],ii[2],field)])/5.0-(wmod[encode3_uin(p,ii[0]-3,ii[1],ii[2],field)]-wmod[encode3_uin(p,ii[0]+3,ii[1],ii[2],field)])/15.0)/2.0))/(2.0*(wd[fencode3_uin(p,ii,delx1)]))    );

  #else
if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[encode3_uin(p,ii[0]+1,ii[1],ii[2],field)]-8*wmod[encode3_uin(p,ii[0]-1,ii[1],ii[2],field)]+wmod[encode3_uin(p,ii[0]-2,ii[1],ii[2],field)]-wmod[encode3_uin(p,ii[0]+2,ii[1],ii[2],field)])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx1)]))    );
 #endif

#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{

  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
  #endif
}
#else

  #ifdef USE_DORDER3
if(ii[0]>2 && ii[0]<((p->n[0])-3) )
 grad=(  ( ((3*wmod[encode3_uin(p,ii[0]+1,ii[1],0,field)]-3*wmod[encode3_uin(p,ii[0]-1,ii[1],0,field)]+3.0*(wmod[encode3_uin(p,ii[0]-2,ii[1],0,field)]-wmod[encode3_uin(p,ii[0]+2,ii[1],0,field)])/5.0-(wmod[encode3_uin(p,ii[0]-3,ii[1],0,field)]-wmod[encode3_uin(p,ii[0]+3,ii[1],0,field)])/15.0)/2.0))/(2.0*(wd[fencode3_uin(p,ii,delx1)]))    );

  #else
if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[encode3_uin(p,ii[0]+1,ii[1],0,field)]-8*wmod[encode3_uin(p,ii[0]-1,ii[1],0,field)]+wmod[encode3_uin(p,ii[0]-2,ii[1],0,field)]-wmod[encode3_uin(p,ii[0]+2,ii[1],0,field)])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx1)]))    );
 #endif
#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{

  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
  #endif
}
#endif



   break;

   case 1:

#ifdef USE_SAC_3D

  #ifdef USE_DORDER3
 if(ii[1]>2 && ii[1]<((p->n[1])-3) )
  grad=(  ( ((3*wmod[encode3_uin(p,ii[0],ii[1]+1,ii[2],field)]-3*wmod[encode3_uin(p,ii[0],ii[1]-1,ii[2],field)]+3.0*(wmod[encode3_uin(p,ii[0],ii[1]-2,ii[2],field)]-wmod[encode3_uin(p,ii[0],ii[1]+2,ii[2],field)])/5.0-(wmod[encode3_uin(p,ii[0],ii[1]-3,ii[2],field)]-wmod[encode3_uin(p,ii[0],ii[1]+3,ii[2],field)])/15.0)/2.0))/(2.0*(wd[fencode3_uin(p,ii,delx2)]))    );

#else
if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[encode3_uin(p,ii[0],ii[1]+1,ii[2],field)]-8*wmod[encode3_uin(p,ii[0],ii[1]-1,ii[2],field)]+wmod[encode3_uin(p,ii[0],ii[1]-2,ii[2],field)]-wmod[encode3_uin(p,ii[0],ii[1]+2,ii[2],field)])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx2)]))    );
 #endif
#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{
  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
  #endif
}
#else

  #ifdef USE_DORDER3
if(ii[1]>2 && ii[1]<((p->n[1])-3) )
 grad=(  ( ((3*wmod[encode3_uin(p,ii[0],ii[1]+1,0,field)]-3*wmod[encode3_uin(p,ii[0],ii[1]-1,0,field)]+3.0*(wmod[encode3_uin(p,ii[0],ii[1]-2,0,field)]-wmod[encode3_uin(p,ii[0],ii[1]+2,0,field)])/5.0-(wmod[encode3_uin(p,ii[0],ii[1]-3,0,field)]-wmod[encode3_uin(p,ii[0],ii[1]+3,0,field)])/15.0)/2.0))/(2.0*(wd[fencode3_uin(p,ii,delx2)]))    );

#endif
if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[encode3_uin(p,ii[0],ii[1]+1,0,field)]-8*wmod[encode3_uin(p,ii[0],ii[1]-1,0,field)]+wmod[encode3_uin(p,ii[0],ii[1]-2,0,field)]-wmod[encode3_uin(p,ii[0],ii[1]+2,0,field)])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx2)]))    );

#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{

  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
  #endif
}
#endif
   break;


   case 2:

#ifdef USE_SAC_3D
  #ifdef USE_DORDER3
 if(ii[2]>2 && ii[2]<((p->n[2])-3) )
  grad=(  ( ((3*wmod[encode3_uin(p,ii[0],ii[1],ii[2]+1,field)]-3*wmod[encode3_uin(p,ii[0],ii[1],ii[2]-1,field)]+3.0*(wmod[encode3_uin(p,ii[0],ii[1],ii[2]-2,field)]-wmod[encode3_uin(p,ii[0],ii[1],ii[2]+2,field)])/5.0-(wmod[encode3_uin(p,ii[0],ii[1],ii[2]-3,field)]-wmod[encode3_uin(p,ii[0],ii[1],ii[2]+3,field)])/15.0)/2.0))/(2.0*(wd[fencode3_uin(p,ii,delx3)]))    );

#else
if( ii[2] >1 &&  ii[2]<((p->n[2])-2))
	grad=(  ( ((8*wmod[encode3_uin(p,ii[0],ii[1],ii[2]+1,field)]-8*wmod[encode3_uin(p,ii[0],ii[1],ii[2]-1,field)]+wmod[encode3_uin(p,ii[0],ii[1],ii[2]-2,field)]-wmod[encode3_uin(p,ii[0],ii[1],ii[2]+2,field)])/6.0))/(2.0*(wd[fencode3_uin(p,ii,delx3)]))    );
#endif

#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{

  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[2]==(p->n[2])-3) || (ii[2]==(p->n[2])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[2]==2 || ii[2]==3  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
  #endif
}
#endif
   break;

}



 return grad;


}









//sum=inte(dpdz,i,j,p->dx[0]);

real inte(real **w,int n, int i, int j, real dx)
{

	real res=0.0;

	if (n == 2)
	  res=dx*0.5*(w[0][j]+w[1][j]);
	else if (n>2)
	{
	  if(i==0) i++;
	  for( int ii=i;ii<n; ii++)
	      res=res+0.5*(w[ii-1][j]+w[ii][j])*dx;


	}

	return res;
}


//bach3d

void initialisation_user1(real *w, real *wd, struct params *p) {




}

void initialisation_user2(real *w, real *wd, struct params *p) {




}





