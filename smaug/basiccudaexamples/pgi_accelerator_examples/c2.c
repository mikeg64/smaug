#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <math.h>
#include <accel.h>
#include <accelmath.h>

int main( int argc, char* argv[] )
{
    int n;      /* size of the vector */
    float *restrict a;  /* the vector */
    float *restrict r;  /* the results */
    float *restrict e;  /* expected results */
    float s, c;
    struct timeval t1, t2, t3;
    long cgpu, chost;
    int i;
    if( argc > 1 )
        n = atoi( argv[1] );
    else
        n = 100000;
    if( n <= 0 ) n = 100000;

    a = (float*)malloc(n*sizeof(float));
    r = (float*)malloc(n*sizeof(float));
    e = (float*)malloc(n*sizeof(float));
    for( i = 0; i < n; ++i ) a[i] = (float)(i+1) * 2.0f;
    /*acc_init( acc_device_nvidia );*/

    gettimeofday( &t1, NULL );
    #pragma acc region
    {
        for( i = 0; i < n; ++i ){
	    s = sinf(a[i]);
	    c = cosf(a[i]);
	    r[i] = s*s + c*c;
	}
    }
    gettimeofday( &t2, NULL );
    cgpu = (t2.tv_sec - t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);
        for( i = 0; i < n; ++i ){
	    s = sinf(a[i]);
	    c = cosf(a[i]);
	    e[i] = s*s + c*c;
	}
    gettimeofday( &t3, NULL );
    chost = (t3.tv_sec - t2.tv_sec)*1000000 + (t3.tv_usec - t2.tv_usec);
    /* check the results */
    for( i = 0; i < n; ++i )
        assert( fabsf(r[i] - e[i]) < 0.000001f );
    printf( "%13d iterations completed\n", n );
    printf( "%13ld microseconds on GPU\n", cgpu );
    printf( "%13ld microseconds on host\n", chost );
    return 0;
}
