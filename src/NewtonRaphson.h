#ifndef _newton
#define _newton

#define MAX_IT 20
#define ACC 0.001

/*
 *
 */
float newton(float x0, float (* f)(float, void *), float (* fdev)(float, void *), void * p);

/*
 * Estimates parameter x of function f in the range [x1,x2]  by means of running the newton iterative algorithm.
 */
float Newton(float x1, float x2, float (* f)(float, void *), float (* fdev)(float, void *), void * p);

#endif
