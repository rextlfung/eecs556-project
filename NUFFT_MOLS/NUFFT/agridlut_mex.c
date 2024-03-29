#include "mex.h" 
#include "ANUFFT2D.c"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])


{

double *sk_real;		
double *sk_imag;		
int K;		

double *sg_real;	
double *sg_imag;	
int N;
int J;	
double *fn;
double *kx1;
double *ky1;
int Ofactor;
int len;
double *dcf;
sk_real = mxGetPr(prhs[0]);		
sk_imag= mxGetPi(prhs[0]);	

K=(int)(*mxGetPr(prhs[1]));


N= (int)(*mxGetPr(prhs[2]));	
            

J= (int)(*mxGetPr(prhs[3]));



fn= mxGetPr(prhs[4]);	
kx1= mxGetPr(prhs[5]);
ky1=mxGetPi(prhs[5]);
Ofactor= (int)(*mxGetPr(prhs[6]));
len=(int)(*mxGetPr(prhs[7]));
dcf=mxGetPr(prhs[8]);
          
plhs[0] = mxCreateDoubleMatrix(1,K*K,mxCOMPLEX);
   
sg_real = mxGetPr(plhs[0]);
sg_imag = mxGetPi(plhs[0]);


gridlut(sk_real,sk_imag,K,N,J,kx1,ky1,sg_real,sg_imag,fn,Ofactor,len,dcf);
		



}
