//comment out one of the init_user routines


//alftest
/*__device__ __host__
void init_user_MODID (real *w, struct params *p,int *ii) {
  int seg1,seg2,seg3,seg4;
  int width=10;

  real m2max=0.001;
  real start=((p->n[0])-width)/2;
  int i,j,k;
  i=ii[0];
  j=ii[1];
  //seg1=((p->n[0])/3)-1;
  seg1=(p->n[0])/6;
  seg2=((p->n[0])/3);
  seg3=(2*(p->n[0])/3)-1;
  //seg4=(2*(p->n[0])/3);
  seg4=(p->n[0])-seg1;



	#ifdef USE_SAC

		    w[fencode3_i(p,ii,rhob)]=1.0;
		    w[fencode3_i(p,ii,b1b)]=1.0;
		    w[fencode3_i(p,ii,energy)]=0.01;

       #else
		    w[fencode3_i(p,ii,rho)]=1.0;
		    w[fencode3_i(p,ii,b1)]=1.0;
		    w[fencode3_i(p,ii,energy)]=0.01;

       #endif

		    //w[fencode3_i(p,ii,b1)]=15*j;
		    //w[fencode3_i(p,ii,b3)]=150*j;
		    
		   //if (i > seg2)
		    //if (i < seg3)
                   // if (i < seg1)
		   //   w[fencode3_i(p,ii,mom2)]=0.0;

		   if (i > seg1)
		    if (i < seg2)
		      w[fencode3_i(p,ii,mom2)]=m2max*(i-seg1)/(seg2-seg1);

		   if (i > seg2)
		    if (i < seg3)
		      //w[fencode3_i(p,ii,mom2)]=m2max*(i-seg2)/(seg3-seg2);
                      w[fencode3_i(p,ii,mom2)]=m2max;
		   if (i > seg3)
		    if (i < seg4)
		      w[fencode3_i(p,ii,mom2)]=m2max*(seg4-i)/(seg4-seg3);






}*/

//ozttest
__device__ __host__
void init_user_MODID (real *w, struct params *p,int *ii) {
                    
                    //real b0=1.0/sqrt(4.0*PI);
                    real b0=1.0;
                    //real ptot=5.0/(12.0*PI);
                    real ptot=5.0/3.0;
                    real rrho=25.0/36.0/PI;
                    real rgamm1;
             real e1,e2;
  int i,j,k;
  i=ii[0];
  j=ii[1];

	#ifdef USE_SAC
                    b0=1.0/sqrt(4.0*PI);
                    ptot=5.0/12.0/PI;
		    w[fencode3_i(p,ii,rhob)]=25.0/(36.0*PI);
                    //w[fencode3_i(p,ii,rhob)]=25.0/9.0;

                    rgamm1=1.0/((p->gamma)-1);
		    
		    w[fencode3_i(p,ii,b1)]=b0*sin((4.0*PI*p->dx[1])*j);
		    w[fencode3_i(p,ii,b2)]=-b0*sin(2.0*PI*(p->dx[0])*i);
		    //w[fencode3_i(p,ii,b2)]=-b0*sin((1.0*p->dx[1])*j);
		    //w[fencode3_i(p,ii,b1)]=b0*sin(2.0*(p->dx[0])*i);

                    //vx=-sin(2pi y)
                    //vy=sin(2pi x)
		    //w[fencode3_i(p,ii,mom1)]=-w[fencode3_i(p,ii,rhob)]*sin(2.0*PI*j*(p->dx[1]));
                    //w[fencode3_i(p,ii,mom2)]=w[fencode3_i(p,ii,rhob)]*sin(2.0*PI*j*(p->dx[0]));

		    w[fencode3_i(p,ii,mom2)]=-w[fencode3_i(p,ii,rhob)]*sin(2.0*PI*i*(p->dx[0]));
                    w[fencode3_i(p,ii,mom1)]=w[fencode3_i(p,ii,rhob)]*sin(2.0*PI*j*(p->dx[1]));
		    //w[fencode3_i(p,ii,mom1)]=-w[fencode3_i(p,ii,rho)]*sin(1.0*i*(p->dx[1]));
                    //w[fencode3_i(p,ii,mom2)]=w[fencode3_i(p,ii,rho)]*sin(1.0*j*(p->dx[0]));

//gives agreement with vac ozt
//w[fencode3_i(p,ii,energy)]=-sin(1.0*i*(p->dx[0]));



                    //w[fencode3_i(p,ii,energy)]=ptot+(0.5*((p->gamma)-1)*(w[fencode3_i(p,ii,mom1)]*w[fencode3_i(p,ii,mom1)]+w[fencode3_i(p,ii,mom2)]*w[fencode3_i(p,ii,mom2)])/rrho);
e1=ptot*rgamm1+(0.5*(w[fencode3_i(p,ii,mom1)]*w[fencode3_i(p,ii,mom1)]+w[fencode3_i(p,ii,mom2)]*w[fencode3_i(p,ii,mom2)])/rrho);
                    //w[fencode3_i(p,ii,energy)]+=0.5*((p->gamma)-2)*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]);

                   e2=0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)])/rrho;
                    //w[fencode3_i(p,ii,energy)]*=rgamm1;
                    w[fencode3_i(p,ii,energyb)]=(e1+e2);

                   w[fencode3_i(p,ii,energy)]=w[fencode3_i(p,ii,energyb)];
                   w[fencode3_i(p,ii,energyb)]=0.0;

                   w[fencode3_i(p,ii,rho)]=w[fencode3_i(p,ii,rhob)];
                   w[fencode3_i(p,ii,rhob)]=0.0;



       #else
		    //w[fencode3_i(p,ii,rho)]=25.0/(36.0*PI);
                    w[fencode3_i(p,ii,rho)]=25.0/9.0;
		    //w[fencode3_i(p,ii,b1)]=-b0*sin((p->dx[0])*i);
		    //w[fencode3_i(p,ii,b2)]=b0*sin(2.0*(p->dx[1])*j);
		    //w[fencode3_i(p,ii,b1)]=b0*sin((2.0*p->dx[1])*(p->n[0] -j));
		    //w[fencode3_i(p,ii,b2)]=-b0*sin(1.0*(p->dx[0])*(p->n[1] -i));
		    w[fencode3_i(p,ii,b1)]=b0*sin((2.0*p->dx[1])*j);
		    w[fencode3_i(p,ii,b2)]=-b0*sin(1.0*(p->dx[0])*i);

		    //w[fencode3_i(p,ii,b3)]=0.0;

                    //vx=-sin(2pi y)
                    //vy=sin(2pi x)
		    //w[fencode3_i(p,ii,mom1)]=-w[fencode3_i(p,ii,rho)]*sin(2.0*PI*j*(p->dx[1]));
                    //w[fencode3_i(p,ii,mom2)]=w[fencode3_i(p,ii,rho)]*sin(2.0*PI*j*(p->dx[0]));
		    //w[fencode3_i(p,ii,mom1)]=-w[fencode3_i(p,ii,rho)]*sin(i*(p->dx[0]));
                    //w[fencode3_i(p,ii,mom2)]=w[fencode3_i(p,ii,rho)]*sin(j*(p->dx[1]));
		    w[fencode3_i(p,ii,mom1)]=w[fencode3_i(p,ii,rho)]*sin(1.0*j*(p->dx[1]));
                    w[fencode3_i(p,ii,mom2)]=-w[fencode3_i(p,ii,rho)]*sin(1.0*i*(p->dx[0]));

		    //w[fencode3_i(p,ii,mom3)]=0;

                    //p=5/12pi  use this to determine the energy
                    //p=(gamma -1)*(e-0.5 rho v**2 - b**2/2)
                    rrho=1.0/w[fencode3_i(p,ii,rho)];
                    rgamm1=1.0/((p->gamma)-1);
		    //w[fencode3_i(p,ii,energy)]=(ptot/((p->gamma)-1))+0.5*rrho*(w[fencode3_i(p,ii,mom1)]*w[fencode3_i(p,ii,mom1)]+w[fencode3_i(p,ii,mom2)]*w[fencode3_i(p,ii,mom2)])+0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]);
//w[fencode3_i(p,ii,energy)]=(ptot/((p->gamma)-1))+0.5*rrho;
		    //w[fencode3_i(p,ii,energy)]=(ptot-(p->gamma)*0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]))*rgamm1+0.5*rrho*(w[fencode3_i(p,ii,mom1)]*w[fencode3_i(p,ii,mom1)]+w[fencode3_i(p,ii,mom2)]*w[fencode3_i(p,ii,mom2)]);
//w[fencode3_i(p,ii,energy)]=(ptot/((p->gamma)-1))+0.5*rrho*(w[fencode3_i(p,ii,mom1)]*w[fencode3_i(p,ii,mom1)]+w[fencode3_i(p,ii,mom2)]*w[fencode3_i(p,ii,mom2)])+0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]);

//gives agreement with vac ozt
w[fencode3_i(p,ii,energy)]=-sin(1.0*i*(p->dx[0]));

       #endif





}

//bwtests
/*__device__ __host__
void init_user_MODID (real *w, struct params *p,int *ii) {
                    
             real p1=1.0;
             real p2=0.1;
             real rho1=1.0;
             real rho2=0.125;
             real rrho=25.0/(36.0*PI);
             real rgamm1;
             real e1,e2;
             int ni=p->n[0];
             int nj=p->n[1];
  int i,j,k;
  i=ii[0];
  j=ii[1];
             e2=0;

	#ifdef USE_SAC
		    

                    rgamm1=1.0/((p->gamma)-1);
		    
                    if(i<(ni*0.315))
                    {

                    w[fencode3_i(p,ii,rhob)]=0.0;
                    w[fencode3_i(p,ii,rho)]=rho1;

		    w[fencode3_i(p,ii,b1)]=0.75;
		    w[fencode3_i(p,ii,b2)]=1.0;
		    w[fencode3_i(p,ii,mom2)]=0.0;
                    w[fencode3_i(p,ii,mom1)]=0.0;
                    //ptot=1.0;
                    //e2=0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]);
                    //e1=(ptot)*rgamm1+(0.5*(w[fencode3_i(p,ii,mom1)]*w[fencode3_i(p,ii,mom1)]+w[fencode3_i(p,ii,mom2)]*w[fencode3_i(p,ii,mom2)])/rrho);
                    e1=(p1)*rgamm1+(0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]));

                    w[fencode3_i(p,ii,energyb)]=(e1);
                    w[fencode3_i(p,ii,energy)]=w[fencode3_i(p,ii,energyb)];
                    w[fencode3_i(p,ii,energyb)]=0.0;


                    }
                    else
                    {

                    w[fencode3_i(p,ii,rhob)]=0.0;
                    w[fencode3_i(p,ii,rho)]=rho2;
 
		    w[fencode3_i(p,ii,b1)]=0.75;
		    w[fencode3_i(p,ii,b2)]=-1.0;
		    w[fencode3_i(p,ii,mom2)]=0.0;
                    w[fencode3_i(p,ii,mom1)]=0.0;
                    //ptot=0.1;
                    //e1=ptot*rgamm1+(0.5*(w[fencode3_i(p,ii,mom1)]*w[fencode3_i(p,ii,mom1)]+w[fencode3_i(p,ii,mom2)]*w[fencode3_i(p,ii,mom2)])/rrho);
                    //e2=0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]);
                    e1=(p2)*rgamm1+(0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]));

                    w[fencode3_i(p,ii,energyb)]=(e1+e2);
                    w[fencode3_i(p,ii,energy)]=w[fencode3_i(p,ii,energyb)];
                    w[fencode3_i(p,ii,energyb)]=0.0;

 

                    }





       #endif





}*/


//bach3d
/*__device__ __host__
void init_user_MODID (real *w, struct params *p,int *ii) {
                    
	real p1,p2,rho0,rho2,v1,v2,v3,T1,T2, xc,yc,zc,r0;
	real Ly, e0,c0;
        real x,y,z;


	Ly=9.46d15;

	e0=1.e48;

	c0=8.95e13;

	p1=1.e0;
	rho0=2.e-22;

	v1=0.e0;
	v2=0.e0;
	v3=0.e0;

	xc=0.0e0;
	yc=0.0e0;
	zc=0.0e0;

	  int i,j,k;
	  i=ii[0];
	  j=ii[1];
	  k=ii[2];

          x=i*(p->dx[0]);
          y=j*(p->dx[1]);
          z=k*(p->dx[2]);
	#ifdef USE_SAC_3D



                    w[fencode3_i(p,ii,rhob)]=0.0;
                    w[fencode3_i(p,ii,energy)]=0.0;
		    w[fencode3_i(p,ii,rhob)]=rho0+c0/((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc));

                    rgamm1=1.0/((p->gamma)-1);
                    w[fencode3_i(p,ii,energyb)]=rgamm1*pow(rho0,(p->gamma));
		    
		    w[fencode3_i(p,ii,b1)]=0;
		    w[fencode3_i(p,ii,b2)]=0;
		    w[fencode3_i(p,ii,b3)]=0;


		    w[fencode3_i(p,ii,mom3)]=v3;
		    w[fencode3_i(p,ii,mom2)]=v2;
		    w[fencode3_i(p,ii,mom1)]=v1;

                   e1=(0.5*rgamm1*(1-((p->gamma)-1))*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]+w[fencode3_i(p,ii,b3)]*w[fencode3_i(p,ii,b3)]));
                    w[fencode3_i(p,ii,energyb)]=w[fencode3_i(p,ii,energyb)]-e1;


                    w[fencode3_i(p,ii,energy)]=w[fencode3_i(p,ii,energyb)];
                    w[fencode3_i(p,ii,energyb)]=0.0;


                     if(i==16 && jj==16  && k==11)
                        w[fencode3_i(p,ii,energy)]=e0/pow(p->dx[0],3.0);
  
			// w(40,28,e_)=e0/(x(1,3,2)-x(1,2,2))**3.d0
			//  w(80,92,e_)=e0/(x(1,3,2)-x(1,2,2))**3.d0  

		    w[fencode3_i(p,ii,bg1)]=w[fencode3_i(p,ii,b1)];
		    w[fencode3_i(p,ii,bg2)]=w[fencode3_i(p,ii,b2)];
		    w[fencode3_i(p,ii,bg3)]=w[fencode3_i(p,ii,b3)];

		    w[fencode3_i(p,ii,b1)]=0;
		    w[fencode3_i(p,ii,b2)]=0;
		    w[fencode3_i(p,ii,b3)]=0;




       #endif





}
*/


