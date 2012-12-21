

__device__ __host__
int addsourceterms2_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

  int direction;
  int status=0;

   real xc1,xc2,xc3;
   real xxmax,yymax;
   real dx,dy,dz;
   real aa;
   real s_period;
   real tdep;

   real vx,vy,vz;

   real exp_x,exp_y,exp_z,exp_xyz;


   real xp,yp,zp;
   int i,j,k;
 	  
	  i=ii[0];
	  j=ii[1];
	  k=ii[2];

     xc1=0.1e6;
    xc2=1.0e6;
    xc3=1.0e6;

          xp=(p->xmin[1])+(((real)j)*(p->dx[1]))-xc2;
          zp=(p->xmin[0])+(((real)i)*(p->dx[0]))-xc1;
          yp=(p->xmin[2])+(((real)k)*(p->dx[2]))-xc3;
     // xx=x(ix_1,ix_2,ix_3,2)-xc2
     // yy=x(ix_1,ix_2,ix_3,3)-xc3
     // zz=x(ix_1,ix_2,ix_3,1)-xc1  


          xp=wd[fencode3_MODID(p,ii,pos2)]-xc2;
          zp=wd[fencode3_MODID(p,ii,pos1)]-xc1;
          yp=wd[fencode3_MODID(p,ii,pos3)]-xc3;  
     

    xxmax=2.0e6;
    yymax=2.0e6;

    dx=0.1e6;
    dy=0.1e6;
    dz=0.05e6;

    aa=10.0;
    s_period=30.0;
    tdep=1.00;


        //exp_z=exp(-zz**2.d0/(delta_z**2.d0))
        //exp_x=exp(-xx**2.d0/(delta_x**2.d0))
        //exp_y=exp(-yy**2.d0/(delta_y**2.d0))       
        //exp_xyz=exp_x*exp_y*exp_z
        exp_z=exp(-zp*zp/(dz*dz));
        exp_x=exp(-xp*xp/(dx*dx));
        exp_y=exp(-yp*yp/(dy*dy));       
        exp_xyz=exp_x*exp_y*exp_z;

        //vvx(ix_1,ix_2,ix_3)=AA*yy/yymax*exp_xyz*tdep    
        //vvy(ix_1,ix_2,ix_3)=-AA*xx/xxmax*exp_xyz*tdep 
        vx=(aa*yp/yymax)*exp_xyz*tdep;    
        vy=-(aa*xp/xxmax)*exp_xyz*tdep;



 

                           w[fencode3_MODID(p,ii,mom2)]+=(p->dt)*vx*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
  
                           w[fencode3_MODID(p,ii,mom3)]+=(p->dt)*vy*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);

                          w[fencode3_MODID(p,ii,energy)]+=(p->dt)*(vx*vx+vy*vy)*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)])/2.0;

        /*if(i==4 && j==63 && k==54) 
	{
                p->test=(w[fencode3_MODID(p,ii,energy)]);
                p->chyp[0]=(w[fencode3_MODID(p,ii,mom2)]);
                p->chyp[1]=(w[fencode3_MODID(p,ii,mom3)]);
                
	}*/

  return ( status);
}

__device__ __host__
int addsourceterms1_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {


}

