
function [data]=NUFFT2_Symmetric(x,K,diff,fn1,kloc,J,N,Ofactor,prefilter_2D,shift)
paddedimage = zeros(K,K);

paddedimage(diff/2+1:end-diff/2,diff/2+1:end-diff/2) = x;

paddedimage = paddedimage.*prefilter_2D;
%KFFTimag=fftshift(fft2(fftshift(paddedimage)));

KFFTimag=fftshift(fft2(fftshift(paddedimage)));
%KFFTimag=fftshift(fft2(paddedimage));


data = giveNUFFT2D(KFFTimag,fn1,kloc,J,K,N,Ofactor);
%data=reshape(data,length(data),1);

if(nargin>9)
    data=data.*shift;
end

end