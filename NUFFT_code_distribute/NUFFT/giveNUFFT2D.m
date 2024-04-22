function data = giveNUFFT2D(KFFTimag,fn1,kloc,J,K,N,Ofactor)



KFFTimag_C=reshape(KFFTimag,1,K*K);


fn = [fn1;0;0;0;0;0];

len=length(kloc);
%%%   kloc is (-N/2,N/2-1)
% zoom K/N

%kloc1=(kloc-1-1i)*(K)/(N);
kloc1=(kloc+N/2+N/2*1i)*(K)/(N);
%kloc1=kloc*K/N;
data=gridlut_mex(KFFTimag_C,K,N,J,fn,kloc1,Ofactor,len);
%data=agridlut_mex(KFFTimag_C,K,N,J,fn,kloc1,Ofactor,len);

end