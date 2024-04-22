function data = giveNUFFT2D_fm(KFFTimag,fn1,fn2,kloc,J,K,N,Ofactor)



KFFTimag_C=reshape(KFFTimag,1,K*K);


fn1 = [0;0;0;0;0;fn1;0;0;0;0;0];
fn2 = [0;0;0;0;0;fn2;0;0;0;0;0];


len=length(kloc);

%kloc 1 to N

kloc1=(kloc-1-1i)*(K)/(N);
%kloc1=kloc*K/N;
data=gridlut2D_mex(KFFTimag_C,K,N,J,fn1,fn2,kloc1,Ofactor,len);

end