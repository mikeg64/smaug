#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "iohyperdif.h"
#include "../include/iotypes.h"

//#include "initialisation_user.h"

#ifdef USE_IOME
void createsim(params k, meta metadata,char *simname, iome el);
void readsim(params *k, meta *md,char *simfile, iome el);
#endif
int encode3_in(Params *dp,int ix, int iy, int iz, int field);
void initconfig(Params *k, Meta *md, real *w, real *wd);
int encode3_uin(Params *dp,int ix, int iy, int iz, int field);
int fencode3_uin (Params *dp,int *ii, int field);
#ifdef USE_SAC_3D
real grad3dngen_uin(real ***wmod, real *wd,struct params *p,int *ii,int dir);
#else
real grad3dngen_uin(real **wmod, real *wd,Params *p,int *ii,int dir);
#endif
real grad3dn_uin(real *wmod, real *wd,struct params *p,int *ii,int field,int dir);
real inte(real **w,int n, int i, int j, real dx);
void initialisation_user1(real *w, real *wd, struct params *p);
void initialisation_user2(real *w, real *wd, struct params *p);
