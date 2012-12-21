#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int encode3_uin(params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
  #else
    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
  #endif
}

int fencode3_uin (struct params *dp,int *ii, int field) {


#ifdef USE_SAC_3D
   return (ii[2]*((dp)->n[0])*((dp)->n[1])  + ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
#else
   return ( ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])));
#endif

}


#ifdef USE_SAC_3D
real grad3dngen_uin(real ***wmod, real *wd,struct params *p,int *ii,int dir)
#else
real grad3dngen_uin(real **wmod, real *wd,struct params *p,int *ii,int dir)
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






