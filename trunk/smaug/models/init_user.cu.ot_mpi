//bwtests
__device__ __host__ void init_user_MODID (real *w, struct params *p,int *ii) {
                    
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

                    rgamm1=1.0/((p->gamma)-1);
                    //swapping dimensions must use the following commented line to
                    //obtain the correct configuration see bug id 540
		    //if((nj-j)<nj*0.315)
                    if(i<(ni*0.315))
                    {
		            w[fencode3_i(p,ii,rhob)]=0.0;
		            w[fencode3_i(p,ii,rho)]=rho1;
			    w[fencode3_i(p,ii,b1)]=0.75;
			    w[fencode3_i(p,ii,b2)]=1.0;
			    w[fencode3_i(p,ii,mom2)]=0.0;
		            w[fencode3_i(p,ii,mom1)]=0.0;
		            //ptot=1.0;
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
		            e1=(p2)*rgamm1+(0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]));
		            w[fencode3_i(p,ii,energyb)]=(e1+e2);
		            w[fencode3_i(p,ii,energy)]=w[fencode3_i(p,ii,energyb)];
		            w[fencode3_i(p,ii,energyb)]=0.0;
                    }
}



