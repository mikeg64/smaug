
//#include "iohyperdif.h"


#include "initialisation_user.h"

#ifdef USE_IOME
void createsim(params k, meta metadata,char *simname, iome el);
void readsim(params *k, meta *md,char *simfile, iome el);
#endif
int encode3_in(params *dp,int ix, int iy, int iz, int field);
void initconfig(params *k, meta *md, real *w, real *wd);



