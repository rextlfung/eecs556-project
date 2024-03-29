

% compute inverse NUFFT 


function gdata = agiveNUFFT2D_fm(data,fn1,fn2,kloc,J,K,N,Ofactor,dcf)


%kloc1=(kloc)*K/N; 
%kloc1 is from 0 to K-1
kloc1=(kloc-1-1i)*(K)/(N);

len=length(kloc1);
fn11 = [0;0;0;0;0;fn1;0;0;0;0;0];
fn12 = [0;0;0;0;0;fn2;0;0;0;0;0];


gdata=agridlut2D_mex(data,K,N,J,fn11,fn12,kloc1,Ofactor,len,dcf);


end