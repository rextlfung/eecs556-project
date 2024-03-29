function Ax = mtimes(A, x)
% spiral_test/MTIMES   Implement Ax for Polynomial&Euclidian decomposition 
% first then b0 field then fft for 3-D data x
% x = [xpeak xpol] where xpeak is the baseline-removed x for only frange; 
% size(mask)=size(x) xpol is the coeffs of x decomposed in polynomial 
% space along f (or n3) dimension
% adapted from fft4b0T2/mtimes


if ~A.adjoint % A*x
    paddedimage = zeros(A.K,A.K);
    diff=A.K-A.N;
    paddedimage(diff/2+1:end-diff/2,diff/2+1:end-diff/2) = x;
    paddedimage = paddedimage.*A.prefilter_2D;    
    KFFTimag=fftshift(fft2(fftshift(paddedimage)))/K;
    data = A.S*KFFTimag(:);
else 
    gdata = reshape(A.St*x,A.K,A.K);
    gdata = fftshift(ifft2(fftshift(gdata)))*K*K/N/N*K;
    gdata=gdata.*prefilter_2D;
    diff=K-N;
    Ax =gdata(diff/2+1:end-diff/2,diff/2+1:end-diff/2); 
end
    
