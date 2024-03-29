
%% Computation of the optimal NUFFT interpolator and prefilter for NUFFT
%
%  The function is NOT assumed to be symmetric and real
%  Computational complexity of this implementation is higher
%  since the kernel is not symmetric or real
%
%  The function computes the optimal mean-square interpolator, which is 
%
% Ref: Z. Yang and M. Jacob,"Mean square optimal NUFFT approximation for efficient non-Cartesian MRI reconstruction" , 
% JMR, submitted
% 
%-------------------------------------------+
% Author: 
% Zhili Yang @ University of Rochester, Mathews Jacob @ University of Iowa
% Email: zhyang@ece.rochester.edu
%        mathews-jacob@uiowa.edu
%-------------------------------------------+
%  Inputs 
%  J: length of the interpolator
%  N: Size of the image
%  K: Oversampling factor of NUFFT. K/N is the oversampling factor
%  Ofactor: Oversampling factor for the lookuptable to store the
%  interpolator; choose >101
%  Order: Order of Bspline interpolation for the lookuptable computation of
%  interpolator
%  Histogram: The energy distribution of the signal; set to ones if unknown
%  

function [prefilter2D,fn,Kernel,error1]=NUFFT_Assymetric(J,N,Ofactor,K,Order,Histogram)

    [fn,Kernel,error1]=  giveAssymetric_MOLS(J,K,N,Ofactor,Order,Histogram);
    prefilter = givePrefilterGeneral(fn,J,K,Ofactor,Order);
    prefilter2D=prefilter*transpose(prefilter);
    
    % scaling the prefilter
   kloc1 = 0+0*i;
   scale =  NUFFT2D_general(ones(N),K,fn,fn,kloc1,J,N,Ofactor,prefilter2D);
   prefilter2D = prefilter2D./scale*N*N;

end