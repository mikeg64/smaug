/////////////////////////////////////
// global variables and configuration section
/////////////////////////////////////

// problem size (vector length) N
static int N = 123456;

// number of threads per block
static int numThreadsPerBlock = 512;

// device to use in case there is more than one
static int selectedDevice = 0;
static int blocksize = 512;

static int dimblock = 16;
