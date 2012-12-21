


//bach3d
__device__ __host__
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
          y=i*(p->dx[1]);
          z=i*(p->dx[2]);
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



