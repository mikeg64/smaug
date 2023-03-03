#include "iotypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int createlog(char *logfile);
int appendlog(char *logfile, Params p, State s);
int writeconfig(char *name,int n,Params p, Meta md, real *w);
int writevtkconfig(char *name,int n,Params p, Meta md, real *w);
int writevacconfig(char *name,int n,Params p, Meta md, real *w,real *wd, State st);
int writevacgatherconfig(char *name,int n,Params p, Meta md, real *w,real *wd, State st);
int readconfig(char *cfgfile, Params p, Meta md, real *w);
int readasciivacconfig(char *cfgfile, Params p, Meta md, State *st, real *w,real *wd, char **hlines, int mode);
/*Big problems with reading fortran unformatted "binary files" need to include
  record field*/
int readbinvacconfig(char *name,Params p, Meta md, real *w,real *wd, State st);
int writeasciivacconfig(char *cfgfile, Params p, Meta md, real *w,real *wd, char **hlines, State st, int mode);
int createconfigsegment(Params p,  real *wnew,real *wdnew, real *w,real *wd);
int gathersegment(Params p,  real *wnew,real *wdnew, real *w,real *wd);
void readatmos(Params p,real *w);
