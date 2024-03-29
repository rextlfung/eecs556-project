function A = NUFFT(K,diff,fn1,kloc,J,N,Ofactor,prefilter_2D,dcf,f)
   A.K=K;
   A.diff=diff;
   A.fn1=fn1;
   A.kloc=kloc;
   A.J=J;
   A.N=N;
   A.Ofactor=Ofactor;
   A.prefilter_2D=prefilter_2D;
   A.adjoint = 0;
   A.dcf=dcf;
   A.f=f;
   A = class(A,'NUFFT');

end
