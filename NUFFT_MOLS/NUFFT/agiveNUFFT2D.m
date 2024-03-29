

% compute inverse NUFFT 


function gdata = agiveNUFFT2D(data,fn1,kloc,J,K,N,Ofactor,dcf)


%kloc1=(kloc)*K/N; 
%kloc1 is from 0 to K-1
kloc1=(kloc+N/2+N/2*1i)*(K)/(N);
%kloc1=(kloc-1-1i)*(K)/(N);
len=length(kloc1);
fn = [fn1;0;0;0;0;0];


gdata=agridlut_mex(data,K,N,J,fn,kloc1,Ofactor,len,dcf);


end