/*  
F = lininterp1f(X,Y,XI,Ydefault) returns the value of the 1-D function Y at the
    points XI using linear interpolation. Length(F)=length(XI).
    The vector X specifies the coordinates of the underlying interval.
    Ydefault is returned for values of XI outside the coordinates in X.
    For lininterp1f to work properly:
    X        must be a monotonically increasing array;
    Y        must be an array with length(Y)=length(X);
	XI       must be an array.
    Ydefault must be a scalar value or an empty matrix [].

	Warning:  not much in the way of error checking, since this slows
              things down, so pay attention to the argument passed to the function!!!

    function lininterp1f(), V. 1.1 Umberto Picchini (umberto.picchini@uniroma1.it), November 2005


    Installation:
	--------------

    Simply copy the trapzf.dll file in a directory recognizable by MATLAB. You won't
	need lininterp1f.c to perform computations, but it is useful if you want to
	read the user instructions or to customize your code.
	If you modify lininterp1f.c then you must run the following command

    mex lininterp1f.c

    in order to obtain the new dll file.


Ex1: 
>> x = [1:1:1000];
>> y =log(sqrt(x+1.001)-1.001);
>> xv =[5:.001:100];
>> yinterp =lininterp1f(x,y,xv,[]);

Ex2: 
>> x=[1:1:1000];
>> y=log(sqrt(x+1.001)-1.001);
>> xv=[5:.001:100];
>> tic; y1=interp1(x,y,xv,'linear'); eval1=toc;  % the MATLAB standard linear interpolator
>> eval1

eval1 =

    0.14

>> tic; y2=lininterp1f(x,y,xv,[]); eval2=toc;
>> eval2

eval2 =

    0.05


Ex3: (run the following script)
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
tic;
rand('seed',0)
for(i=1:50)
    x = [1:1:1000] + rand(1,1000);
    y = x.^(1./x);
    xv = [3:.001:100];
    yv = interp1(x,y,xv,'linear');   % the MATLAB standard linear interpolator
end
eval1 = toc;

tic;
rand('seed',0)
for(i=1:50)
    x = [1:1:1000] + rand(1,1000);
    y = x.^(1./x);
    xv = [3:.001:100];
    yv = lininterp1f(x,y,xv,[]);
end
eval2 = toc;
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


then you obtain:

eval1 =

     9.1640


eval2 =

    3.4750
*/


#include "mex.h"

/* Input arguments */
#define X_data      prhs[0]
#define Y_data      prhs[1]
#define X_interp    prhs[2]
#define Y_default   prhs[3]

/* Output arguments */
#define Y_out       plhs[0]



static double lininterp1f(double *yinterp, double *xv, double *yv, double *x, double *ydefault, int m, int minterp)
{
    int i, j; 
	int nrowsinterp, nrowsdata;
	nrowsinterp = minterp;
	nrowsdata = m;
	for (i=0; i<nrowsinterp; i++)
	{
			if((x[i] < xv[0]) || (x[i] > xv[nrowsdata-1]))
				yinterp[i] = *ydefault;
			else
			{   for(j=1; j<nrowsdata; j++)
			{      if(x[i]<=xv[j])
			{		   yinterp[i] = (x[i]-xv[j-1]) / (xv[j]-xv[j-1]) * (yv[j]-yv[j-1]) + yv[j-1];
                       break;
			}
			}
			}
	}
    return *yinterp;
} 


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int m, minterp;
double *yinterp;
double *xv, *yv, *x, *ydefault;

/* Find the dimension of the data */
   m = mxGetNumberOfElements(X_data);  
/* Find the dimension of the x values for which interpolation is desired */
   minterp = mxGetNumberOfElements(X_interp);
/* Create an mxArray of real values for the output */
   Y_out = mxCreateDoubleMatrix(1,minterp,mxREAL);
/* Get the data passed in */
   xv       = mxGetPr(X_data);
   yv       = mxGetPr(Y_data);
   x        = mxGetPr(X_interp);
   ydefault = mxGetPr(Y_default);
/* the output */
   yinterp  = mxGetPr(Y_out);
/* Do the actual computations in a subroutine */
   lininterp1f(yinterp,xv,yv,x,ydefault,m,minterp);
   return;
}