#include "mex.h" 
#include "NUFFT2D_fm.c"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])


{

 
	
double *sg_real;		
double *sg_imag;		
int K;		

double *sk_real;	
double *sk_imag;	
int N;
int J;	
double *fnx_real;
double *fnx_imag;
double *fny_real;
double *fny_imag;
double *kx1;
double *ky1;

int Ofactor;
int len;

sg_real = mxGetPr(prhs[0]);		
sg_imag= mxGetPi(prhs[0]);	

if(sg_imag==NULL)
  mexPrintf("no imaginary part!"); 

K=(int)(*mxGetPr(prhs[1]));


N= (int)(*mxGetPr(prhs[2]));	
            

J= (int)(*mxGetPr(prhs[3]));



fnx_real= mxGetPr(prhs[4]);	
fnx_imag= mxGetPi(prhs[4]);
fny_real= mxGetPr(prhs[5]);	
fny_imag= mxGetPi(prhs[5]);

kx1= mxGetPr(prhs[6]);
ky1= mxGetPi(prhs[6]);

Ofactor= (int)(*mxGetPr(prhs[7]));
len=(int)(*mxGetPr(prhs[8]));
          
plhs[0] = mxCreateDoubleMatrix(1,len,mxCOMPLEX);
   
sk_real = mxGetPr(plhs[0]);
sk_imag = mxGetPi(plhs[0]);


gridlut2d(sg_real,sg_imag,K,N,J,kx1,ky1,sk_real,sk_imag,fnx_real,fnx_imag,fny_real,fny_imag,Ofactor,len);
		



}
