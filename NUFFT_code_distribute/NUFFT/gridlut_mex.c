#include "mex.h" 
#include "NUFFT2.c"
//#include "ANUFFT2D.c"
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])


{

 
	
double *sg_real;		
double *sg_imag;		
int K;		

double *sk_real;	
double *sk_imag;	
int N;
int J;	
double *fn;
double *kx1;
double *ky1;
int Ofactor;
int len;

sg_real = mxGetPr(prhs[0]);		
sg_imag= mxGetPi(prhs[0]);	

//if(s_imag==NULL)
  // mexPrintf("no imaginary part!"); 

K=(int)(*mxGetPr(prhs[1]));


N= (int)(*mxGetPr(prhs[2]));	
            

J= (int)(*mxGetPr(prhs[3]));



fn= mxGetPr(prhs[4]);	
kx1= mxGetPr(prhs[5]);
ky1=mxGetPi(prhs[5]);
Ofactor= (int)(*mxGetPr(prhs[6]));
len=(int)(*mxGetPr(prhs[7]));
          
plhs[0] = mxCreateDoubleMatrix(1,len,mxCOMPLEX);
   
sk_real = mxGetPr(plhs[0]);
sk_imag = mxGetPi(plhs[0]);


gridlut(sg_real,sg_imag,K,N,J,kx1,ky1,sk_real,sk_imag,fn,Ofactor,len);
		



}

