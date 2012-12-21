#include "iotypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int createlog(char *logfile);
int appendlog(char *logfile, params p, state s);
int writeconfig(char *name,int n,params p, meta md, real *w);
int writevtkconfig(char *name,int n,params p, meta md, real *w);
int writevacconfig(char *name,int n,params p, meta md, real *w,real *wd, state st);
int writevacgatherconfig(char *name,int n,params p, meta md, real *w,real *wd, state st);
int readconfig(char *cfgfile, params p, meta md, real *w);
int readasciivacconfig(char *cfgfile, params p, meta md, state *st, real *w,real *wd, char **hlines, int mode);
/*Big problems with reading fortran unformatted "binary files" need to include 
  record field*/
int readbinvacconfig(char *name,params p, meta md, real *w,real *wd, state st);
int writeasciivacconfig(char *cfgfile, params p, meta md, real *w,real *wd, char **hlines, state st, int mode);
int createconfigsegment(params p,  real *wnew,real *wdnew, real *w,real *wd);
int gathersegment(params p,  real *wnew,real *wdnew, real *w,real *wd);
void readatmos(params p,real *w);
