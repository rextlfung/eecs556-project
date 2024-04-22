function idata=INUFFT2D_general(data,fn1,fn2,kloc1,J,K,N,Ofactor,dcf,prefilter2D_fm,shift)
    
kloc=kloc1+N/2+N/2*1i+1+1i;
%kloc = kloc1;

if(nargin>10)
    data = data.*conj(shift);
end
    gdata = agiveNUFFT2D_fm(data,fn1,fn2,kloc,J,K,N,Ofactor,dcf);
    gdata=reshape(gdata,K,K);
    gdata=fftshift(ifft2(fftshift(gdata)));
    gdata=gdata.*prefilter2D_fm;
    diff=K-N;
    idata=gdata(diff/2+1:end-diff/2,diff/2+1:end-diff/2); 

end