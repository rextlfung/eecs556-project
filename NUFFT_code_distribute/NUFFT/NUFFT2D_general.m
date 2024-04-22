% %
%-------------------------------------------+
% Author: 
% Zhili Yang @ University of Rochester, Mathews Jacob @ University of Iowa
% Email: zhyang@ece.rochester.edu
%        mathews-jacob@uiowa.edu
%-------------------------------------------+
 
%% Description of the parameters
%
% x – Original image, scaled to [0,1];
% N - Size of image, here is 128;
% K - Oversampling factor, classical option is K=2N. here due to memory
% concern, we prefer K/N much smaller than 2, close to 1.
% J - Size of support limited [-J/2,J/2] interpolators. 
% fn1 and fn2 are interpolators along x and y directions
% kloc - Non-cartesian trajectory
% X - Simulated k-space data
% Ofactor- interpolator oversampling factor. 
% prefilter_2D -scale factors

%%


function [data]=NUFFT2D_general(x,K,fn1,fn2,kloc1,J,N,Ofactor,prefilter_2D,shift)

    kloc=kloc1+N/2+N/2*1i+1+1i;
%    kloc=kloc1;

    paddedimage = zeros(K,K);
    diff=K-N;
    paddedimage(diff/2+1:end-diff/2,diff/2+1:end-diff/2) = x;

    paddedimage = paddedimage.*prefilter_2D;

    KFFTimag=fftshift(fft2(fftshift(paddedimage)));
    data = giveNUFFT2D_fm(KFFTimag,fn1,fn2,kloc,J,K,N,Ofactor);
    
    if(nargin>9)
        data=data.*shift;
    end

end