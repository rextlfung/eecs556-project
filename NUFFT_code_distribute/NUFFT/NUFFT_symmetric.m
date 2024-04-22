%-------------------------------------------+
% Author: 
% Zhili Yang @ University of Rochester, Mathews Jacob @ University of Iowa
% Email: zhyang@ece.rochester.edu
%        mathews-jacob@uiowa.edu
%-------------------------------------------+
 

%% Function to compute the optimal least square symmetric interpolator and the associated prefilter
%
%  The function is assumed to be symmetric and real
%  Computational complexity of this implementation is low
%
% Ref: M. Jacob,"Optimized least square non uniform fast Fourier transform (OLS-NUFFT)" , 
% IEEE Transactions of Signal Processing, vol. 57, issue 6, pp. 2165-2177, Feb 2009 
%
% J - Size of support limited [-J/2,J/2] interpolators. 
% N - Size of image, here is 128;
% K - Oversampling factor, classical option is K=2N. 
% Ofactor- interpolator oversampling factor, used for lookup table based storage of interpolator. 
% Order - Order of Bspline used for lookup table interpolation

%%
function [prefilter_2D,fn,fnfull,kernel,error]=NUFFT_symmetric(J,N,Ofactor,K,Order)


    % compute symmetric worst case optimal interpolating function
    [fn,error,kernel]=  giveSymmetricNUFFTfn(J,K,N,Ofactor,Order);  
    
    fnfull = [flipud(fn(2:end));fn];
    
    % compute prefilter/scale factors
    prefilter = givePrefilter(fnfull,J,K,Ofactor,Order);
    prefilter_2D=prefilter*prefilter';
   
    % scaling the prefilter
    kloc = 0+0*i;
    scale = NUFFT2_Symmetric(ones(N),K,K-N,fn,kloc,J,N,Ofactor,prefilter_2D);
    prefilter_2D = prefilter_2D./scale*N*N;

end