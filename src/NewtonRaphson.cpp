#include "NewtonRaphson.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>


using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;


float newton(float x0, float (* f)(float, void *), float (* fdev)(float, void *), void * p)
{
    int i;
    float xi;
    float xi1;
    float dx;

    xi1 = xi = x0;
    //cout<<"x0: "<<x0<<endl;
    
    float fvalue,fdevalue;
    for (i = 0; i < MAX_IT; i++)
    {
        fvalue=f(xi, p);
        fdevalue=fdev(xi, p);
        xi1 = xi - (fvalue/fdevalue);
        
      //  cout<<"xi1: "<<xi1<<" f: "<<fvalue<<" fdev: "<<fdevalue<<endl;
        dx = fabs(xi1 - xi);
        xi = xi1;
        if (dx < ACC)
        {
            break;
        }
    }

    return xi;
}

/*
 * Estimates parameter x of function f in the range [x1,x2]  by means of running the newton iterative algorithm.
 */
float Newton(float x1, float x2, float (* f)(float, void *), float (* fdev)(float, void *), void * p)
{
  //  cout<<"x1: "<<x1<<" x2: "<<x2<<endl;
    int j;
    float df,dx,dxold,ff,fh,fl;
    float temp,xh,xl,rts;

    fl = f(x1, p);
    df = fdev(x1, p);
    fh = f(x2, p);
    df = fdev(x2, p);

    if (fl == 0.0)
        return x1;
    if (fh == 0.0)
        return x2;
    if (fl < 0.0)
    {
        xl=x1;
        xh=x2;
    }
    else
    {
        xh=x1;
        xl=x2;
    }
    rts=0.5*(x1+x2);
    rts=x1;//<-fijamos un valor inicial para el x.
    dxold=fabs(x2-x1);
    dx=dxold;

    ff = f(rts, p);
    df = fdev(rts, p);
   
    for (j=1;j<=MAX_IT;j++)
    {
     //   cout<<"rts: "<<rts<<endl;
        //Initialize the guess for root, the “stepsize before last,” and the last step.
        //Loop over allowed iterations.
        if ((((rts-xh)*df-ff)*((rts-xl)*df-ff) > 0.0) || (fabs(2.0*ff) > fabs(dxold*df)))	//Bisect if Newton out of range,
        {	//or not decreasing fast enough.
   
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            if (xl == rts)
                return rts;
        }
        else
        {
   
            dxold=dx;
            dx=ff/df;
            temp=rts;
            rts -= dx;
            if (temp == rts)
                return rts;
            //Change in root is negligible. Newton step acceptable. Take it.
        }
        if (fabs(dx) < ACC){
            return rts;
        }
        ff = f(rts, p);
        df = fdev(rts, p);
        //Orient the search so that f(xl) < 0.
        //Convergence criterion.

        if (ff < 0.0)	//Maintain the bracket on the root.
            xl=rts;
        else
            xh=rts;
    }

    return 10.0;
}
