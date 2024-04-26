
function [idata, Xk] =INUFFT2_Symmetric(data,fn1,kloc,J,K,N,Ofactor,dcf,prefilter_2D,shift);


if(nargin>9)
    data = data.*conj(shift);
end

gdata = agiveNUFFT2D(data,fn1,kloc,J,K,N,Ofactor,dcf);
gdata=reshape(gdata,K,K);

Xk = gdata;

gdata=fftshift(ifft2(fftshift(gdata)));

gdata=gdata.*prefilter_2D;

diff=K-N;
idata=gdata(diff/2+1:end-diff/2,diff/2+1:end-diff/2); 

end