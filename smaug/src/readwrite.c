#include "../include/readwrite.h"

void freadl(FILE *stream, char **string)
{    
    unsigned long counter = 0;
    char *line = NULL;
    int next = fgetc(stream);

    do {
        next = fgetc(stream);
        if (next == EOF) {
            free(line);
            break;
        }
        ++counter;
        line = (char*)realloc(line, counter + 1);
        if (line == NULL) {
            puts("line == NULL");
            exit(EXIT_FAILURE);
        }
        line[counter - 1] = (char)next;
    } while (next != '\n');
    line[counter - 1] = '\0';
    *string = (char *)malloc(strlen(line) + 1);
    if (*string == NULL) {
        puts("*string == NULL");
    } else {
        strcpy(*string, line);
    }
    free(line);

}



int createlog(char *logfile)
{
	int status=0;

      FILE *fdt=0;

      fdt=fopen(logfile,"a+");
      fprintf(fdt,"it   t   dt    rho m1 m2 m3 e bx by bz\n");
      fclose(fdt);	

	return status;
}

int appendlog(char *logfile, params p, state s)
{
  int status=0;
      FILE *fdt=0;

      fdt=fopen(logfile,"a+");
      fprintf(fdt,"%d %f %f %f %f %f %f %f %f %f %f %f\n",s.it,
               s.t,s.dt,s.rho,s.m1,s.m2,s.m3,s.e,s.b1,s.b2,s.b3);
      fclose(fdt);
  return status;
}

int writeconfig(char *name,int n,params p, meta md, real *w)
{
  int status=0;
  int i1,j1;
  int ni,nj;
  char configfile[300];


  ni=p.n[0];
  nj=p.n[1];



      //save file containing current data
      sprintf(configfile,"tmp/%ss%d.out",name,n);
      printf("check dims %d %d \n",ni,nj);
      FILE *fdt=fopen(configfile,"w");
      fprintf(fdt,"%d\n",n);
     for( j1=(p.ng[1]);j1<(nj-(p.ng[1]));j1++)
      {
        for( i1=(p.ng[0]);i1<(nj-(p.ng[0]));i1++)
	{
               // printf("%d %d ", i1,j1);
 #ifdef ADIABHYDRO
		fprintf(fdt,"%d %d %f %f %f %f \n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*energy)]);
#endif
 #ifdef USE_SAC
		fprintf(fdt,"%d %d %f %f %f %f %f %f\n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)]);
#endif
 #ifdef USE_SAC_3D
		fprintf(fdt,"%d %d %f %f %f %f %f %f %f %f\n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*mom3)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)],w[j1*ni+i1+(ni*nj*b3)]);
#endif
           //fprintf(fdt,"%d %f %f %f ",j1+i1*nj, u[j1+i1*nj],v[j1+i1*nj],h[j1+i1*nj]);
               // fprintf(fdt,"%f ",h[j1+i1*nj]);
        }     
        //printf("\n");   
        //fprintf(fdt,"\n");
      }
      fclose(fdt);


      //save file containing current data
      sprintf(configfile,"out/%s.out",name);
      printf("write out check dims %s %d %d \n",configfile ,ni,nj);
      fdt=fopen(configfile,"a+");
      fprintf(fdt,"%d\n",n);
     for( j1=(p.ng[1]);j1<(nj-(p.ng[1]));j1++)
      {
        for( i1=(p.ng[0]);i1<(nj-(p.ng[0]));i1++)
	{
               // printf("%d %d ", i1,j1);
 #ifdef ADIABHYDRO
 fprintf(fdt,"%d %d %f %f %f %f \n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*energy)]);
#endif
#ifdef USE_SAC
fprintf(fdt,"%d %d %f %f %f %f %f %f\n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)]);
        #endif
#ifdef USE_SAC_3D
fprintf(fdt,"%d %d %f %f %f %f %f %f %f %f\n",i1,j1,w[(j1*ni+i1)+(ni*nj*rho)],w[(j1*ni+i1)+(ni*nj*mom1)],w[(j1*ni+i1)+(ni*nj*mom2)],w[j1*ni+i1+(ni*nj*mom3)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)],w[j1*ni+i1+(ni*nj*b3)]);
        #endif
		
           //fprintf(fdt,"%d %f %f %f ",j1+i1*nj, u[j1+i1*nj],v[j1+i1*nj],h[j1+i1*nj]);
               // fprintf(fdt,"%f ",h[j1+i1*nj]);
        }     
        //printf("\n");   
        //fprintf(fdt,"\n");
      }
      fclose(fdt);


  return status;
}



int writevacconfig(char *name,int n,params p, meta md, real *w, real *wd, state st)
{
  int status=0;
  int i1,j1,k1,ifield;
  int ni,nj,nk;
  
  char tcfg[300];  
  char configfile[300];
  
  char buffer[800];
  double dbuffer[12];
  int ibuffer[5];
  char ext[3];

  ni=p.n[0];
  nj=p.n[1];
    #ifdef USE_SAC_3D
  nk=p.n[2];
    #endif

 
   sprintf(configfile,"%s",name);
   
   #ifdef USE_MPI
   
   char *pch1,*pch2;
   pch1 = strtok (configfile,".");
   sprintf(tcfg,"%s",pch1);
   pch2 = strtok (NULL,".");

   printf("here1 %s %d %d \n",tcfg,ni,nj);
   //sprintf(ext,"%s",pch2);

      //set the input filename corresponding to proc id
 
     #ifdef USE_SAC3D
	      if(p.ipe>99)
		sprintf(configfile,"%s%d_np%d%d%d_%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe);
	      else if(p.ipe>9)
		sprintf(configfile,"%s%d_np0%d0%d0%d_0%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe);
	      else
		sprintf(configfile,"%s%d_np00%d00%d00%d_00%d.out",tcfg,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe);  	     
     #else
	      if(p.ipe>99)
		sprintf(configfile,"%s%d_np%d%d_%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.ipe);
	      else if(p.ipe>9)
		sprintf(configfile,"%s%d_np%d%d_%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.ipe);
	      else
		sprintf(configfile,"%s%d_np0%d0%d_00%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.ipe);  	     	     
     #endif

 
   #else


         //save file containing current data
      sprintf(configfile,"%s_%d.out",name,st.it);
  #endif
   

     // sprintf(configfile,"%s",name);
      printf("write vac check dims %d %d %d %lf\n",ni,nj,st.it,st.t);
//printf("here3 %s \n",configfile); 
      FILE *fdt=fopen(configfile,"w");
      //FILE *fdt=fopen("out/test.out","w");


      fwrite(md.name,sizeof(char)*79,1,fdt);

      // ibuffer[0]=24;
       //fwrite(ibuffer,sizeof(int),1,fdt);
      //fwrite(md.name,sizeof(char)*79,1,fdt);
      //*line2:
      //*   it          - timestep (integer)
      //*   t           - time     (real)
      //*   ndim        - dimensionality, negative sign for gen. coord (integer)
      //*   neqpar      - number of equation parameters (integer)
      //*   nw          - number of flow variables (integer)
      //sprintf(buffer,"%ld %lg %ld %ld %ld\n",st.it,st.t,3,4,8);
      //it,time,ndim,neqpar,nw
      printf("st.it=%d\n",st.it);
      ibuffer[0]=st.it;
      dbuffer[0]=st.t;
      //ibuffer[0]=1;
      //dbuffer[0]=10.98; st.it
      fwrite(ibuffer,sizeof(int),1,fdt);
      fwrite(dbuffer,sizeof(double),1,fdt);

    #ifdef USE_SAC
      ibuffer[0]=2;
      ibuffer[1]=6;
      ibuffer[2]=10;
      //ibuffer[0]=st.it;
      //ibuffer[1]=(float)st.t;


      //ibuffer[2]=2;  //report from distribution routine to reset the data
      //ibuffer[3]=6;
      //ibuffer[4]=10;
   #endif
    #ifdef USE_SAC_3D
      ibuffer[0]=3;
      ibuffer[1]=7;
      ibuffer[2]=13;
      //ibuffer[2]=3;
    #endif
      fwrite(ibuffer,sizeof(int)*3,1,fdt);
      //fwrite(ibuffer,sizeof(int)*5,1,fdt);
      //line3:
      //*   nx()        - the grid dimensions      (ndim integers)
      //sprintf(buffer,"%ld %ld\n",ni,nj);
      ibuffer[0]=ni;
      ibuffer[1]=nj;

      //ibuffer[0]=128;
      //ibuffer[1]=256;
    #ifdef USE_SAC
      fwrite(ibuffer,sizeof(int)*2,1,fdt);
    #endif
    #ifdef USE_SAC_3D
      ibuffer[2]=nk;
      fwrite(ibuffer,sizeof(int)*3,1,fdt);
    #endif
      //*line4:
      //*   eqpar()     - equation parameters from filenameini (neqpar reals)
      //sprintf(buffer,"%lg %lg %lg %lg %lg %lg\n",p.gamma,p.eta,p.g[0],p.g[1],0,0);
      dbuffer[0]=p.gamma;
      dbuffer[1]=p.eta;
      dbuffer[2]=p.g[0];
      dbuffer[3]=p.g[1];

    #ifdef USE_SAC
      dbuffer[4]=0;
      dbuffer[5]=0;
     fwrite(dbuffer,sizeof(double)*6,1,fdt);
    #endif
    #ifdef USE_SAC_3D
      dbuffer[4]=p.g[2];
      dbuffer[5]=0;
      dbuffer[6]=0;
     fwrite(dbuffer,sizeof(double)*7,1,fdt);
    #endif
 

      //*line5:
      //*   varnames    - names of the coordinates, variables, equation parameters
      //*                 eg. 'x y rho mx my e bx by  gamma eta' (character*79)

    #ifdef USE_SAC_3D
      sprintf(buffer,"x y z rho mx my mz e bx by bz gamma eta g1 g2 g3");
    #else
      sprintf(buffer,"x y rho mx my e bx by  gamma eta g1 g2");
    #endif

      fwrite(buffer,sizeof(char)*79,1,fdt);

    #ifdef USE_SAC_3D
      for(ifield=0;ifield<16;ifield++)   
   #else
       for(ifield=0;ifield<12;ifield++)   
   #endif 

    #ifdef USE_SAC_3D
   for( k1=0;k1<nk;k1++)
    #endif
for( j1=0;j1<nj;j1++)
  
	{
//energyb,rhob,b1b,b2b         
       for( i1=0;i1<ni;i1++)     
      {
          //if(ifield==2) printf("density %lG ",w[(j1*ni+i1)]);
               if(ifield==0)
               dbuffer[0]=wd[(j1*ni+i1)+(ni*nj*(pos1))];
               else if(ifield==1)
               dbuffer[0]=wd[(j1*ni+i1)+(ni*nj*(pos2))];
    #ifdef USE_SAC_3D
               if(ifield==0)
               dbuffer[0]=wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos1))];
               else if(ifield==1)
               dbuffer[0]=wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos2))];
               else if(ifield==2)
               dbuffer[0]=wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos3))];
               
    #else
    
    #endif
               else
    #ifdef USE_SAC_3D
                dbuffer[0]=w[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(ifield-3))];
    #else
                dbuffer[0]=w[(j1*ni+i1)+(ni*nj*(ifield-2))];              
    #endif

                fwrite(dbuffer,sizeof(double),1,fdt);
                //printf("%g ",dbuffer[0]);		

        }     
      }
      //printf("density nowt \n \n ");
      buffer[0]='\n';
      fwrite(buffer,sizeof(char),1,fdt);
     // printf("\n");
      fclose(fdt);

  return status;
}


int writevacgatherconfig(char *name,int n,params p, meta md, real *w, real *wd, state st)
{
  int status=0;
  int i1,j1,k1,ifield;
  int ni,nj,nk;
  int istart,jstart,kstart;
  int ifin,jfin,kfin;
  
  char tcfg[300];  
  char configfile[300];
  
  char buffer[800];
  double dbuffer[12];
  int ibuffer[5];
  char ext[3];

 istart=2*((p.noghost)>0);
 jstart=2*((p.noghost)>0);
 kstart=2*((p.noghost)>0);



  ni=p.n[0];//-4*((p.noghost)>0);
  nj=p.n[1];//-4*((p.noghost)>0);

 ifin=ni-2*((p.noghost)>0);
 jfin=nj-2*((p.noghost)>0);

    #ifdef USE_SAC_3D
  nk=p.n[2];//-4*((p.noghost)>0);
  kfin=nk+2*((p.noghost)>0);
    #endif




 
   sprintf(configfile,"%s",name);
      sprintf(configfile,"%s",name);
   
   #ifdef USE_MPI
   
   char *pch1,*pch2;
   pch1 = strtok (configfile,".");
   sprintf(tcfg,"%s",pch1);
   pch2 = strtok (NULL,".");

   printf("here1 %s %d %d \n",tcfg,ni,nj);
     sprintf(configfile,"%s%d.out",tcfg,st.it);
   #else
 


         //save file containing current data
      sprintf(configfile,"%s_%d.out",name,st.it);
  #endif
   

     // sprintf(configfile,"%s",name);
      printf("write vac check dims %d %d %d %lf\n",ni,nj,st.it,st.t);
//printf("here3 %s \n",configfile); 
      FILE *fdt=fopen(configfile,"w");
      //FILE *fdt=fopen("out/test.out","w");


      fwrite(md.name,sizeof(char)*79,1,fdt);

      // ibuffer[0]=24;
       //fwrite(ibuffer,sizeof(int),1,fdt);
      //fwrite(md.name,sizeof(char)*79,1,fdt);
      //*line2:
      //*   it          - timestep (integer)
      //*   t           - time     (real)
      //*   ndim        - dimensionality, negative sign for gen. coord (integer)
      //*   neqpar      - number of equation parameters (integer)
      //*   nw          - number of flow variables (integer)
      //sprintf(buffer,"%ld %lg %ld %ld %ld\n",st.it,st.t,3,4,8);
      //it,time,ndim,neqpar,nw
      printf("st.it=%d\n",st.it);
      ibuffer[0]=st.it;
      dbuffer[0]=st.t;
      //ibuffer[0]=1;
      //dbuffer[0]=10.98; st.it
      fwrite(ibuffer,sizeof(int),1,fdt);
      fwrite(dbuffer,sizeof(double),1,fdt);

    #ifdef USE_SAC
      ibuffer[0]=2;
      ibuffer[1]=6;
      ibuffer[2]=10;
      //ibuffer[0]=st.it;
      //ibuffer[1]=(float)st.t;


      //ibuffer[2]=2;  //report from distribution routine to reset the data
      //ibuffer[3]=6;
      //ibuffer[4]=10;
   #endif
    #ifdef USE_SAC_3D
      ibuffer[0]=3;
      ibuffer[1]=7;
      ibuffer[2]=13;
      //ibuffer[2]=3;
    #endif
      fwrite(ibuffer,sizeof(int)*3,1,fdt);
      //fwrite(ibuffer,sizeof(int)*5,1,fdt);
      //line3:
      //*   nx()        - the grid dimensions      (ndim integers)
      //sprintf(buffer,"%ld %ld\n",ni,nj);
      ibuffer[0]=ni-4*((p.noghost)>0);
      ibuffer[1]=nj-4*((p.noghost)>0);

      //ibuffer[0]=128;
      //ibuffer[1]=256;
    #ifdef USE_SAC
      fwrite(ibuffer,sizeof(int)*2,1,fdt);
    #endif
    #ifdef USE_SAC_3D
      ibuffer[2]=nk-4*((p.noghost)>0);
      fwrite(ibuffer,sizeof(int)*3,1,fdt);
    #endif
      //*line4:
      //*   eqpar()     - equation parameters from filenameini (neqpar reals)
      //sprintf(buffer,"%lg %lg %lg %lg %lg %lg\n",p.gamma,p.eta,p.g[0],p.g[1],0,0);
      dbuffer[0]=p.gamma;
      dbuffer[1]=p.eta;
      dbuffer[2]=p.g[0];
      dbuffer[3]=p.g[1];

    #ifdef USE_SAC
      dbuffer[4]=0;
      dbuffer[5]=0;
     fwrite(dbuffer,sizeof(double)*6,1,fdt);
    #endif
    #ifdef USE_SAC_3D
      dbuffer[4]=p.g[2];
      dbuffer[5]=0;
      dbuffer[6]=0;
     fwrite(dbuffer,sizeof(double)*7,1,fdt);
    #endif
 

      //*line5:
      //*   varnames    - names of the coordinates, variables, equation parameters
      //*                 eg. 'x y rho mx my e bx by  gamma eta' (character*79)

    #ifdef USE_SAC_3D
      sprintf(buffer,"x y z rho mx my mz e bx by bz gamma eta g1 g2 g3");
    #else
      sprintf(buffer,"x y rho mx my e bx by  gamma eta g1 g2");
    #endif

      fwrite(buffer,sizeof(char)*79,1,fdt);

    #ifdef USE_SAC_3D
      for(ifield=0;ifield<16;ifield++)   
   #else
       for(ifield=0;ifield<12;ifield++)   
   #endif 

    #ifdef USE_SAC_3D
   for( k1=kstart;k1<kfin;k1++)
    #endif
for( j1=jstart;j1<jfin;j1++)
  
	{
//energyb,rhob,b1b,b2b         
       for( i1=istart;i1<ifin;i1++)     
      {
          //if(ifield==2) printf("density %lG ",w[(j1*ni+i1)]);
               if(ifield==0)
               dbuffer[0]=wd[(j1*ni+i1)+(ni*nj*(pos1))];
               else if(ifield==1)
               dbuffer[0]=wd[(j1*ni+i1)+(ni*nj*(pos2))];
    #ifdef USE_SAC_3D
               if(ifield==0)
               dbuffer[0]=wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos1))];
               else if(ifield==1)
               dbuffer[0]=wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos2))];
               else if(ifield==2)
               dbuffer[0]=wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos3))];
               
    #else
    
    #endif
               else
    #ifdef USE_SAC_3D
                dbuffer[0]=w[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(ifield-3))];
    #else
                dbuffer[0]=w[(j1*ni+i1)+(ni*nj*(ifield-2))];              
    #endif

                fwrite(dbuffer,sizeof(double),1,fdt);		
                //printf("%g ",w[(j1*ni+i1)+(ni*nj*(ifield-2))]);	
        }     
      }
      //printf("density nowt \n \n ");
      buffer[0]='\n';
      //printf("\n ");	
      fwrite(buffer,sizeof(char),1,fdt);
      fclose(fdt);

  return status;
}




/*Big problems with reading fortran unformatted "binary files" need to include 
  record field*/

int readbinvacconfig(char *name,params p, meta md, real *w,real *wd, state st)
{

  int status=0;
  int i1,j1,k1,ifield;
  int ni,nj,nk;
  
  char tcfg[300];  
  char configfile[300];
  
  char buffer[800];
  double dbuffer[12];
  int ibuffer[5];
  char ext[3];
  long lsize;
  size_t result;

  char *bigbuf;
  ni=p.n[0];
  nj=p.n[1];
    #ifdef USE_SAC_3D
  nk=p.n[2];
    #endif






      //save file containing current data
      //sprintf(configfile,"out/v%s.out",name);
      sprintf(configfile,"%s",name);
      printf("check dims %d %d \n",ni,nj);
      FILE *fdt=fopen(configfile,"r");


      fread(buffer,sizeof(char),80,fdt);
      for(i1=0;i1<81;i1++)
         putchar(buffer[i1]);
      printf("starting %s\n ",buffer);
       strcpy(md.name,buffer);

      //*line2:
      //*   it          - timestep (integer)
      //*   t           - time     (real)
      //*   ndim        - dimensionality, negative sign for gen. coord (integer)
      //*   neqpar      - number of equation parameters (integer)
      //*   nw          - number of flow variables (integer)
      //sprintf(buffer,"%ld %lg %ld %ld %ld\n",st.it,st.t,3,4,8);
      //it,time,ndim,neqpar,nw

      fread(ibuffer,sizeof(int),1,fdt);
      fread(dbuffer,sizeof(double),1,fdt);
      st.it=ibuffer[0]=st.it;
      st.t=dbuffer[0]=st.t;
      printf("st.it=%f st.t=%d\n",st.it,st.t);
      fread(ibuffer,sizeof(int),3,fdt);

      //fread(ibuffer,sizeof(int)*3,1,fdt);

      //line3:
      //*   nx()        - the grid dimensions      (ndim integers)
      //sprintf(buffer,"%ld %ld\n",ni,nj);
      //fread(ibuffer,sizeof(int)*2,1,fdt);
    #ifdef USE_SAC
      fread(ibuffer,sizeof(int)*2,1,fdt);
    #endif
    #ifdef USE_SAC_3D
      
      fread(ibuffer,sizeof(int)*3,1,fdt);
      //nk=ibuffer[2];
    #endif
      //ni=ibuffer[0];
      //nj=ibuffer[1];







      //*line4:
      //*   eqpar()     - equation parameters from filenameini (neqpar reals)
      //sprintf(buffer,"%lg %lg %lg %lg %lg %lg\n",p.gamma,p.eta,p.g[0],p.g[1],0,0);
      //fread(dbuffer,sizeof(double)*6,1,fdt);



    #ifdef USE_SAC
      fread(dbuffer,sizeof(double)*6,1,fdt);
      //dbuffer[4]=0;
      //dbuffer[5]=0;
     
    #endif
    #ifdef USE_SAC_3D
      fread(dbuffer,sizeof(double)*7,1,fdt);
      p.g[2]=dbuffer[4];
      //dbuffer[5]=0;
      //dbuffer[6]=0;     
    #endif

      p.gamma=dbuffer[0];
      p.eta=dbuffer[1];
      p.g[0]=dbuffer[2];
      p.g[1]=dbuffer[3];
      printf("%f %f %f %f\n",dbuffer[0],dbuffer[1],dbuffer[2],dbuffer[3]);



      //*line5:
      //*   varnames    - names of the coordinates, variables, equation parameters
      //*                 eg. 'x y rho mx my e bx by  gamma eta' (character*79)
      sprintf(buffer,"x y rho mx my mz e bx by bz gamma eta g1 g2 g3\n");
      fread(buffer,sizeof(char)*79,1,fdt);



      // printf("read buffer%s\n",buffer);



    #ifdef USE_SAC_3D
      for(ifield=0;ifield<16;ifield++)   
   #else
       for(ifield=0;ifield<12;ifield++)   
   #endif 

    #ifdef USE_SAC_3D
   for( k1=0;k1<nk;k1++)
    #endif


         for( i1=0;i1<ni;i1++)   
	{
//energyb,rhob,b1b,b2b         
     
for( j1=0;j1<nj;j1++)  
      {



               fread(dbuffer,sizeof(double),1,fdt);
               if(ifield==0)
               	wd[(j1*ni+i1)+(ni*nj*(pos1))]=dbuffer[0];
               else if(ifield==1)
               	wd[(j1*ni+i1)+(ni*nj*(pos2))]=dbuffer[0];
    #ifdef USE_SAC_3D
               if(ifield==0)
                 wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos1))]=dbuffer[0];
               else if(ifield==1)
                 wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos2))]=dbuffer[0];
               else if(ifield==2)
                 wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos3))]=dbuffer[0];
               
    #else
    
    #endif
               else
    #ifdef USE_SAC_3D
                w[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(ifield-3))]=dbuffer[0];
    #else
                w[(j1*ni+i1)+(ni*nj*(ifield-2))]=dbuffer[0]; 
              //if(ifield==2)
              //  printf("%g ",w[(j1*ni+i1)+(ni*nj*(ifield-2))]);             
    #endif

               






               // fread(dbuffer,12*sizeof(double),1,fdt);		
                //i1*p.dx[0]=dbuffer[0];
                //j1*p.dx[1]=dbuffer[1];
               //if(j1==2 || j1==3)
                //  printf("%d %d %d %d %lg\n", ni,nj,i1,j1,dbuffer[2]);
 


        } 
        // if(ifield==2)
       // printf("\n");    
      }

      // printf("read bin vac read fields\n");
      fclose(fdt);
      free(bigbuf);
  return status;


}


int writevtkconfig(char *name,int n,params p, meta md, real *w)
{
  int status=0;
  int i1,j1,k1;
  int ni,nj,nk;
  char configfile[300];
  char labels[4][4]={"rho","e","mom","b"};
  int is;
  ni=p.n[0];
  nj=p.n[1];

#ifdef USE_SAC_3D

nk=p.n[2];
#endif


      //save file containing current data

      //scalar fields
//n+=10;
        #ifdef ADIABHYDRO
      for(int i=0,is=0; i<=3; i+=4)
      //for(int i=0; i<=4; i+=3)
       #else
      for(int i=0,is=0; i<=4; i+=4,is+=3)
      //for(int i=0; i<=4; i+=3)
        #endif
      {
	      if(n<=9)
                 sprintf(configfile,"%s%ss00%d.vtk",name,labels[i/4],n);
              else if(n<=99)
                 sprintf(configfile,"%s%ss0%d.vtk",name,labels[i/4],n);
              else
                 sprintf(configfile,"%s%ss%d.vtk",name,labels[i/4],n);

	      //printf("check dims %s %s %d %d \n",configfile,labels[i/4],ni,nj);
	      FILE *fdt=fopen(configfile,"w");


	      fprintf(fdt,"# vtk DataFile Version 2.0\n");
	      fprintf(fdt,"Structured Grid\n");
	      fprintf(fdt,"ASCII\n");
	      fprintf(fdt," \n");
	      fprintf(fdt,"DATASET RECTILINEAR_GRID\n");
#ifdef USE_SAC_3D
	      fprintf(fdt,"DIMENSIONS %d %d %d\n",ni,nj,nk);
#else
	      fprintf(fdt,"DIMENSIONS %d %d 1\n",ni,nj);
#endif


	      fprintf(fdt,"X_COORDINATES %d double\n",ni);
              for(i1=0;i1<ni;i1++)
	        fprintf(fdt,"%G\n",(p.xmin[0])+i1*p.dx[0]);

	      fprintf(fdt,"Y_COORDINATES %d double\n",nj);
              for(i1=0;i1<nj;i1++)
	        fprintf(fdt,"%G\n",(p.xmin[1])+i1*p.dx[1]);


   #ifdef USE_SAC_3D
	      fprintf(fdt,"Z_COORDINATES %d double\n",nk);
              for(k1=0;k1<nk;k1++)
	        fprintf(fdt,"%G\n",(p.xmin[2])+k1*p.dx[2]);

    #else
	      fprintf(fdt,"Z_COORDINATES 1 double\n");
	      fprintf(fdt,"0\n");
   #endif
   #ifdef USE_SAC_3D
	      fprintf(fdt,"POINT_DATA  %d\n",(ni)*(nj)*nk);
   #else

	      fprintf(fdt,"POINT_DATA  %d\n",(ni)*(nj));
    #endif

	      fprintf(fdt,"SCALARS %s double 1\n",labels[i/4]);

             fprintf(fdt,"LOOKUP_TABLE TableName \n");

   #ifdef USE_SAC_3D
	     for( k1=0;k1<(nk);k1++)
   #endif
	     for( j1=0;j1<(nj);j1++)
		for( i1=0;i1<(ni);i1++)
                {
                 if(is==0)
                    #ifdef USE_SAC_3D
			fprintf(fdt,"%G\n",w[(k1*ni*nj)+(j1*ni+i1)+(ni*nj*nk*is)]+w[(k1*ni*nj)+(j1*ni+i1)+(ni*nj*nk*8)]);
                    #else
			fprintf(fdt,"%G\n",w[(j1*ni+i1)+(ni*nj*is)]+w[(j1*ni+i1)+(ni*nj*7)]);
                    #endif
                 else
                    #ifdef USE_SAC_3D
			fprintf(fdt,"%G\n",w[(k1*ni*nj)+(j1*ni+i1)+(ni*nj*nk*is)]);
                    #else
			fprintf(fdt,"%G\n",w[(j1*ni+i1)+(ni*nj*is)]);
                    #endif

                }

	      fclose(fdt);
      }

      //vector fields
      int iv;
        #ifdef ADIABHYDRO
      for(int i=2; i<3; i++)
       #else
      for(int i=2; i<=3; i++)
        #endif

      {
	      if(i==2)
                iv=1;
              else
                //iv=5;
                iv=4;
              if(n<=9)
                 sprintf(configfile,"vtk/%s%ss00%d.vtk",labels[i],name,n);
              else if(n<=99)
                 sprintf(configfile,"vtk/%s%ss0%d.vtk",labels[i],name,n);
              else
                 sprintf(configfile,"vtk/%s%ss%d.vtk",labels[i],name,n);

	      //printf("check dims %s %s %d %d \n",configfile,labels[i],ni,nj);
	      FILE *fdt=fopen(configfile,"w");


	      fprintf(fdt,"# vtk DataFile Version 2.0\n");
	      fprintf(fdt,"Structured Grid\n");
	      fprintf(fdt,"ASCII\n");
	      fprintf(fdt," \n");
	      fprintf(fdt,"DATASET RECTILINEAR_GRID\n");
	      fprintf(fdt,"DIMENSIONS %d %d 1\n",ni,nj);


	      fprintf(fdt,"X_COORDINATES %d double\n",ni);
              for(i1=0;i1<ni;i1++)
	        fprintf(fdt,"%G\n",(p.xmin[0])+i1*p.dx[0]);

	      fprintf(fdt,"Y_COORDINATES %d double\n",nj);
              for(i1=0;i1<nj;i1++)
	        fprintf(fdt,"%G\n",(p.xmin[1])+i1*p.dx[1]);


               #ifdef USE_SAC_3D
	      fprintf(fdt,"Z_COORDINATES %d double\n",nk);
              for(i1=0;i1<nk;i1++)
	        fprintf(fdt,"%G\n",(p.xmin[2])+i1*p.dx[2]);
               #else
	      fprintf(fdt,"Z_COORDINATES 1 double\n");
	      fprintf(fdt,"0\n");
               #endif


             #ifdef USE_SAC_3D
	      fprintf(fdt,"POINT_DATA  %d\n",(ni)*(nj)*nk);
             #else
	      fprintf(fdt,"POINT_DATA  %d\n",(ni)*(nj));
             #endif
	      fprintf(fdt,"VECTORS %s double \n",labels[i]);

            #ifdef USE_SAC_3D
		for( k1=0;k1<(nk);k1++)
             #endif
		for( j1=0;j1<(nj);j1++)
	      		for( i1=0;i1<(ni);i1++)
   
            #ifdef USE_SAC_3D
                         fprintf(fdt,"%G %G %G\n",w[(k1*ni*nj)+(j1*ni+i1)+(ni*nk*nj*iv)],w[(k1*ni*nj)+(j1*ni+i1)+(ni*nk*nj*(iv+1))],w[(k1*ni*nj)+(j1*ni+i1)+(ni*nk*nj*(iv+2))]);
            #else
                         fprintf(fdt,"%G %G %G\n",w[(j1*ni+i1)+(ni*nj*iv)],w[(j1*ni+i1)+(ni*nj*(iv+1))]);
             #endif

                       //printing mag fields including backround for SAC
                       //if(iv==4)
                       //  fprintf(fdt,"%f %f %f\n",w[(j1*ni+i1)+(ni*nj*iv)]+w[(j1*ni+i1)+(ni*nj*(iv+4))],w[(j1*ni+i1)+(ni*nj*(iv+1))]+w[(j1*ni+i1)+(ni*nj*(iv+1+4))]);


	      fclose(fdt);
      }




  return status;
}



int readconfig(char *cfgfile, params p, meta md, real *w)
{
  int status=0;

  return status;
}


int readasciivacconfig(char *cfgfile, params p, meta md,state *st, real *w, real *wd, char **hlines, int mode)
{
  int status=0;
  int i;
  int i1,j1,k1;
  int ni,nj,nk;
  int iif,jf,kf;
  int is,js,ks;
  int shift;
  real x,y,z,val;
  char cfgfilename[300];
  char ext[3];

   int ii1,ii2,ii3;

   ni=p.n[0];
   nj=p.n[1];
   is=0;
   js=0;
   ks=0;

   iif=ni;
   jf=nj;

   sprintf(cfgfilename,"%s",cfgfile);
   #ifdef USE_SAC_3D
   nk=p.n[2];
   kf=nk;
   #endif

    if(mode==0)
    {
       #ifdef USE_MULTIGPU
		if((p.pipe[0])==0) iif=ni -2;
		if((p.pipe[0])==((p.pnpe[0])-1)) is=2;
		if((p.pipe[1])==0) jf=nj-2;
		if((p.pipe[1])==((p.pnpe[1])-1)) js=2;

	  #ifdef USE_SAC_3D
		if((p.pipe[2])==0) kf=nk-2;
		if((p.pipe[2])==((p.pnpe[2])-1)) ks=2;
          #endif
	
        #endif
    }
  
printf("reading %s %d %d\n",cfgfile,ni,nj);
   FILE *fdt=fopen(cfgfilename,"r+");
//FILE *fdt=fopen("zero1_np0201_001.ini","r+");
   //char **hlines;
   char *line;
   //hlines=(char **)calloc(5, sizeof(char*));
freadl(fdt, &hlines[0]);
     printf("%s\n", hlines[0]);

fscanf(fdt,"%d %lG %d %d %d\n",&(st->it),&(st->t),&ii1,&ii2,&ii3);

   //read 5 header lines
   for(i=2;i<5;i++)
   {
     freadl(fdt, &hlines[i]);
     printf("%s\n", hlines[i]);
   }
//printf("read header %d\n" , p.ipe);
  //fscanf(fdt,"%f",&val);
 //printf("%f",val);
#ifdef USE_SAC_3D
for( k1=ks;k1<(kf);k1++)
#endif
for( j1=js;j1<(jf);j1++)
for( i1=is;i1<(iif);i1++)

   
	     
             {

#ifdef USE_SAC_3D
                         shift=(k1*ni*nj+j1*ni+i1);
                         fscanf(fdt,"%lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG\n",&wd[shift+(ni*nj*nk*pos1)],&wd[shift+(ni*nj*nk*pos2)],&wd[shift+(ni*nj*nk*pos3)], &w[shift],&w[shift+(ni*nj*nk)],&w[shift+(ni*nj*nk*2)],&w[shift+(ni*nj*nk*3)],&w[shift+(ni*nj*nk*4)],&w[shift+(ni*nj*nk*5)],&w[shift+(ni*nj*nk*6)],&w[shift+(ni*nj*nk*7)],&w[shift+(ni*nj*nk*8)],&w[shift+(ni*nj*nk*9)],&w[shift+(ni*nj*nk*10)],&w[shift+(ni*nj*nk*11)],&w[shift+(ni*nj*nk*12)]);


#else
                         shift=(j1*ni+i1);
                         fscanf(fdt,"%lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG\n",&wd[shift+(ni*nj*pos1)],&wd[shift+(ni*nj*pos2)],&w[shift],&w[shift+(ni*nj)],&w[shift+(ni*nj*2)],&w[shift+(ni*nj*3)],&w[shift+(ni*nj*4)],&w[shift+(ni*nj*5)],&w[shift+(ni*nj*6)],&w[shift+(ni*nj*7)],&w[shift+(ni*nj*8)],&w[shift+(ni*nj*9)]);

//if(p.ipe==0  && w[shift]<=0)
//printf("density %lG %lG %lG \n",wd[shift+(ni*nj*pos1)],wd[shift+(ni*nj*pos1)],w[shift]);


#endif

                         //freadl(fdt, &line);
                         /*for(i=0; i<NVAR;i++)
                         {
                            //fscanf(fdt,"%g",&);
                             //w[shift+(ni*nj*i)]=val;
                              printf("%lG ",w[shift+(ni*nj*i)]);
                         }*/
                         //fscanf(fdt,"\n");
                         //printf("\n");
              }


	      fclose(fdt);

  //free(hlines);
  return status;
}

int writeasciivacconfig(char *cfgfile, params p, meta md, real *w,real *wd, char **hlines, state st, int mode)
{
  int status=0;
  int i;
  int i1,j1;
  int ni,nj;                         
  int shift;
  real x,y,val;

  int iif,jf,kf;
  int is,js,ks;
  is=0;
  js=0;
  ks=0;

   ni=p.n[0];
   nj=p.n[1];

   iif=ni;
   jf=nj;
    
	  #ifdef USE_SAC_3D
            nk=p.n[2];
            kf=nk;
	#endif
   
   //char **hlines;
   char *line;

 char configfile[300];
  char tcfg[300];

sprintf(configfile,"%s",cfgfile);

if(mode==0)
{
      #ifdef USE_MULTIGPU
		if((p.pipe[0])==0) iif=ni -2;
		if((p.pipe[0])==((p.pnpe[0])-1)) is=2;
		if((p.pipe[1])==0) jf=nj-2;
		if((p.pipe[1])==((p.pnpe[1])-1)) js=2;

	  #ifdef USE_SAC_3D
		if((p.pipe[2])==0) kf=nk-2;
		if((p.pipe[2])==((p.pnpe[2])-1)) ks=2;
          #endif
	
        #endif
}


 #ifdef USE_MPI
   char ext[3];
   char *pch1,*pch2;
   pch1 = strtok (configfile,".");
   sprintf(tcfg,"%s",pch1);
   pch2 = strtok (NULL,".");
   sprintf(ext,"%s",pch2);


 

   //printf("here1 %s \n",tcfg);
   //sprintf(ext,"%s",pch2);

      //set the input filename corresponding to proc id
      //if scattering configuration and ini file will not have the name appended with the iteration number
     if(strcmp(ext,"ini")==0)
     {
	     #ifdef USE_SAC3D
		      if(p.ipe>99)
			sprintf(configfile,"%s_np%d%d%d_%d.%s",tcfg,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe,ext);
		      else if(p.ipe>9)
			sprintf(configfile,"%s_np0%d0%d0%d_0%d.%s",tcfg,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe,ext);
		      else
			sprintf(configfile,"%s_np00%d00%d00%d_00%d.%s",tcfg,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe,ext);  	     
	     #else
		      if(p.ipe>99)
			sprintf(configfile,"%s_np%d%d_%d.%s",tcfg,p.pnpe[0],p.pnpe[1],p.ipe,ext);
		      else if(p.ipe>9)
			sprintf(configfile,"%s_np%d%d_%d.%s",tcfg,p.pnpe[0],p.pnpe[1],p.ipe,ext);
		      else
			sprintf(configfile,"%s_np0%d0%d_00%d.%s",tcfg,p.pnpe[0],p.pnpe[1],p.ipe,ext);  	     	     
	     #endif

     }
     else
     {
	     #ifdef USE_SAC3D
		      if(p.ipe>99)
			sprintf(configfile,"%s%d_np%d%d%d_%d.%s",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe,ext);
		      else if(p.ipe>9)
			sprintf(configfile,"%s%d_np0%d0%d0%d_0%d.%s",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe,ext);
		      else
			sprintf(configfile,"%s%d_np00%d00%d00%d_00%d.%s",tcfg,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe,ext);  	     
	     #else
		      if(p.ipe>99)
			sprintf(configfile,"%s%d_np%d%d_%d.%s",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.ipe,ext);
		      else if(p.ipe>9)
			sprintf(configfile,"%s%d_np%d%d_%d.%s",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.ipe,ext);
		      else
			sprintf(configfile,"%s%d_np0%d0%d_00%d.%s",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.ipe,ext);  	     	     
	     #endif
     }

    // if(p.mode==3)
    //    sprintf(configfile,"%s",cfgfile);
   #else
   char ext[3];
   char *pch1,*pch2;
   pch1 = strtok (configfile,".");
   sprintf(tcfg,"%s",pch1);
   pch2 = strtok (NULL,".");
   sprintf(ext,"%s",pch2);
     if(strcmp(ext,"ini")==0)
     {
         sprintf(configfile,"%s.%s",tcfg,ext);
       }
     else
     {
	sprintf(configfile,"%s_%d.%s",tcfg,st.it,ext); 
	//sprintf(configfile,"%s_%d.out",cfgfile,st.it);
       }

         //save file containing current data
      
  #endif

//printf("%s %d", configfile,p.ipe);
//return 0;
  FILE *fdt=fopen(configfile,"a+");



   //hlines=(char **)calloc(5, sizeof(char*));

   //printf("here %s\n",hlines[0]);
   //fprintf(fdt,"%s\n", hlines[0]);
   fprintf(fdt,"%s\n", md.name);
   //read 5 header lines

      //*line2:
      //*   it          - timestep (integer)
      //*   t           - time     (real)
      //*   ndim        - dimensionality, negative sign for gen. coord (integer)
      //*   neqpar      - number of equation parameters (integer)
      //*   nw          - number of flow variables (integer)
     // fprintf(fdt,"%d %d %d 6 %d\n", st.it,st.t,NDIM,NVAR);

    #ifdef USE_SAC_3D
     fprintf(fdt,"%d %lg %d 7 %d\n", st.it,st.t,NDIM,NVAR);
     printf("%d %g %d 7 %d\n", st.it,st.t,NDIM,NVAR);
    #else
     fprintf(fdt,"%d %lg %d 6 %d\n", st.it,st.t,NDIM,NVAR);
     printf("ascii write %d %g %d 6 %d\n", st.it,st.t,NDIM,NVAR);
   #endif

       #ifdef USE_SAC
      fprintf(fdt,"%d %d\n",ni,nj);
    #endif
    #ifdef USE_SAC_3D
      fprintf(fdt,"%d %d %d\n",ni,nj,nk);
    #endif 

     //*line4:
      //*   eqpar()     - equation parameters from filenameini (neqpar reals)
      //sprintf(buffer,"%lg %lg %lg %lg %lg %lg\n",p.gamma,p.eta,p.g[0],p.g[1],0,0);

    #ifdef USE_SAC
      fprintf(fdt,"%lg %lg %lg %lg %d %d\n",p.gamma,p.eta,p.g[0],p.g[1],0,0);
    #endif
    #ifdef USE_SAC_3D
      fprintf(fdt,"%lg %lg %lg %lg %lg %d %d\n",p.gamma,p.eta,p.g[0],p.g[1],p.g[2],0,0);
    #endif

      //*line5:
      //*   varnames    - names of the coordinates, variables, equation parameters
      //*                 eg. 'x y rho mx my e bx by  gamma eta' (character*79)

    #ifdef USE_SAC_3D
      fprintf(fdt,"x y z rho mx my mz e bx by bz gamma eta g1 g2 g3\n");
    #else
      fprintf(fdt,"x y rho mx my e bx by  gamma eta g1 g2\n");
    #endif


#ifdef USE_SAC_3D
for( k1=ks;k1<(kf);k1++)
#endif
   for( j1=js;j1<(jf);j1++)
	     for( i1=is;i1<(iif);i1++)
             {
                         x=(1+i1)*(p.dx[0]);
                         y=(1+j1)*(p.dx[1]);
			#ifdef USE_SAC_3D
				z=(1+k1)*(p.dx[2]);
			#endif

                         //shift=(j1*ni+i1);
                         //fprintf(fdt,"%lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE\n",wd[shift+ni*nj*pos1],wd[shift+ni*nj*pos2],w[shift],w[shift+(ni*nj)],w[shift+(ni*nj*2)],w[shift+(ni*nj*3)],w[shift+(ni*nj*4)],w[shift+(ni*nj*5)],w[shift+(ni*nj*6)],w[shift+(ni*nj*7)],w[shift+(ni*nj*8)],w[shift+(ni*nj*9)]);

                         //printf(" %g %d ",w[shift],shift);

#ifdef USE_SAC_3D
                         shift=(k1*ni*nj+j1*ni+i1);
                         fprintf(fdt,"%lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG\n",wd[shift+(ni*nj*nk*pos1)],wd[shift+(ni*nj*nk*pos2)],wd[shift+(ni*nj*nk*pos3)], w[shift],w[shift+(ni*nj*nk)],w[shift+(ni*nj*nk*2)],w[shift+(ni*nj*nk*3)],w[shift+(ni*nj*nk*4)],w[shift+(ni*nj*nk*5)],w[shift+(ni*nj*nk*6)],w[shift+(ni*nj*nk*7)],w[shift+(ni*nj*nk*8)],w[shift+(ni*nj*nk*9)],w[shift+(ni*nj*nk*10)],w[shift+(ni*nj*nk*11)],w[shift+(ni*nj*nk*12)]);


#else
                         shift=(j1*ni+i1);
                         fprintf(fdt,"%lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG %lG\n",wd[shift+(ni*nj*pos1)],wd[shift+(ni*nj*pos2)],w[shift],w[shift+(ni*nj)],w[shift+(ni*nj*2)],w[shift+(ni*nj*3)],w[shift+(ni*nj*4)],w[shift+(ni*nj*5)],w[shift+(ni*nj*6)],w[shift+(ni*nj*7)],w[shift+(ni*nj*8)],w[shift+(ni*nj*9)]);

//if(p.ipe==0)
//printf("density %lG %lG %lG \n",x,y,w[shift]);



#endif






              }


	      fclose(fdt);

  //free(hlines);
  return status;
}


int createconfigsegment(params p,  real *wnew,real *wdnew, real *w,real *wd)
{
  int status=0;
  int i,var;
  int i1,j1,k1;
  int ni,nj;                         
  int shift,tvar;
  real x,y,val;

  int oni,onj,onk;

   k1=0;

   #ifdef USE_MULTIGPU
   ni=p.n[0]/(p.pnpe[0]);
   nj=p.n[1]/(p.pnpe[1]);
   oni=p.n[0];
   onj=p.n[1];

int oi1,oj1,ok1;
int oshift;

   #ifdef USE_SAC3D
   nk=p.n[2]/(p.pnpe[2]);
   onk=p.n[2];


//for(tvar=0; tvar<NVAR+NDIM; tvar++)
   for( k1=0;k1<(nk);k1++)

   #endif
   for( j1=0;j1<(nj);j1++)
   {
	     for( i1=0;i1<(ni);i1++)
             {
                oi1=i1+(p.pipe[0]*ni);
                oj1=j1+(p.pipe[1]*nj);  
		#ifdef USE_SAC_3D
                         shift=(k1*ni*nj+j1*ni+i1);
                         ok1=k1+(p.pipe[2]*nk);

                         oshift=(ok1*oni*onj+oj1*oni+oi1);
		#else
			 shift=(j1*ni+i1);
                         oshift=(oj1*oni+oi1);
                #endif
                   for(var=0; var<NVAR; var++)
                    {
                         wnew[shift+ni*nj*var]=w[oshift+oni*onj*var];

                       /*  if(var ==2)
                            printf("%g ",w[oshift+oni*onj*var]);*/
                     }
                 

                    for(var=pos1; var<=(pos2+(NDIM>2)); var++)
                         wdnew[shift+ni*nj*var]=wd[oshift+oni*onj*var];
//                         fprintf(fdt,"%lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE\n",wd[shift+ni*nj*pos1],wd[shift+ni*nj*pos2],w[shift],w[shift+(ni*nj)],w[shift+(ni*nj*2)],w[shift+(ni*nj*3)],w[shift+(ni*nj*4)],w[shift+(ni*nj*5)],w[shift+(ni*nj*6)],w[shift+(ni*nj*7)],w[shift+(ni*nj*8)],w[shift+(ni*nj*9)]);




              }
       // printf("\n");

   }

  #endif
  //free(hlines);
  return status;
}


int gathersegment(params p,  real *wnew,real *wdnew, real *w,real *wd)
{
  int status=0;
  int i,var;
  int i1,j1,k1;
  int ni,nj;                         
  int shift;
  real x,y,val;

  int oni,onj,onk;


   #ifdef USE_MULTIGPU

   k1=0;
   ni=p.n[0];
   nj=p.n[1];
   oni=ni*(p.pnpe[0]);
   onj=nj*(p.pnpe[1]);

int oi1,oj1,ok1;
int oshift;

   #ifdef USE_SAC3D
   nk=p.n[2];
   onk=nk*(p.pnpe[2]);

   for( k1=0;k1<(nk);k1++)

   #endif
   for( j1=0;j1<(nj);j1++)
	     for( i1=0;i1<(ni);i1++)
             {
                oi1=i1+(p.pipe[0]*ni);
                oj1=j1+(p.pipe[1]*nj);  
		#ifdef USE_SAC_3D
                         shift=(k1*ni*nj+j1*ni+i1);
                         ok1=k1+(p.pipe[2]*nk);

                         oshift=(ok1*oni*onj+oj1*oni+oi1);
		#else
			 shift=(j1*ni+i1);
                         oshift=(oj1*oni+oi1);
                #endif
                    for(var=0; var<NVAR; var++)
                         wnew[oshift+oni*onj*var]=w[shift+ni*nj*var];

                    for(var=pos1; var<=(pos2+(NDIM>2)); var++)
                         wdnew[oshift+oni*onj*var]=wd[shift+ni*nj*var];
//                         fprintf(fdt,"%lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE\n",wd[shift+ni*nj*pos1],wd[shift+ni*nj*pos2],w[shift],w[shift+(ni*nj)],w[shift+(ni*nj*2)],w[shift+(ni*nj*3)],w[shift+(ni*nj*4)],w[shift+(ni*nj*5)],w[shift+(ni*nj*6)],w[shift+(ni*nj*7)],w[shift+(ni*nj*8)],w[shift+(ni*nj*9)]);




              }

  #endif
  //free(hlines);
  return status;
}


void readatmos(params p,real *w)
{
int n=-1;
int count;
real h,rho0;
char s1[100], s2[100];
FILE *fatmos;
if(fatmos=fopen("test.dat","r"))
{
//fprintf(fatmos,"%f %f %n\n",1.3,1.4);
//rewind(fatmos);

for(int i=0;i<4;i++)
{
count=fscanf(fatmos, " %s %s %n", s1, s2,&n);
//freadl_ui(fatmos,&inbuf);
if(-1 == n) printf("parse failed\n");
printf("fatmos read   %s %s %d\n",s1,s2,count);
}
fclose(fatmos);
}
else
{
printf("failed\n");
}





}




