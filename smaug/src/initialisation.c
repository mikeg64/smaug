
#include <stdio.h>
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
