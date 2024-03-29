function Ax = mtimes(A, x)
% R_gradxy/MTIMES   Implement Rx for computing xy-gradient of 4-D data "x"
% and also its adjoint

if ~A.adjoint % R*x
    Ax = NUFFT2_v(x,A.K,A.diff,A.fn1,A.kloc,A.J,A.N,A.Ofactor,A.prefilter_2D,A.f); % 
else % R'*x
    Ax = INUFFT2_v(x,A.fn1,A.kloc,A.J,A.K,A.N,A.Ofactor,A.dcf,A.prefilter_2D,A.f); %  
end
end

function [data1]=NUFFT2_v(x,K,diff,fn1,kloc,J,N,Ofactor,prefilter_2D,f)
paddedimage = zeros(K,K);

paddedimage(diff/2+1:end-diff/2,diff/2+1:end-diff/2) = x;

paddedimage = paddedimage.*prefilter_2D;

KFFTimag=fftshift(fft2(fftshift(paddedimage)));
%KFFTimag=fft2(paddedimage);
data = giveNUFFT2D(KFFTimag,fn1,kloc,J,K,N,Ofactor);
%data1=data;

    data1=data.*f;

end

function [idata]=INUFFT2_v(data,fn1,kloc,J,K,N,Ofactor,dcf,prefilter_2D,f);

data1=data./f;
gdata = agiveNUFFT2D(data1,fn1,kloc,J,K,N,Ofactor,dcf);



gdata=reshape(gdata,K,K);


gdata=fftshift(ifft2(fftshift(gdata)));
%gdata=fftshift(ifft2(gdata));

gdata=gdata.*prefilter_2D;


diff=K-N;
idata=gdata(diff/2+1:end-diff/2,diff/2+1:end-diff/2); 


end
% 
%     function [data1]=NUFFT2_v(x,K,diff,fn1,kloc,J,N,Ofactor,prefilter_2D)
% paddedimage = zeros(K,K);
% 
% paddedimage(diff/2+1:end-diff/2,diff/2+1:end-diff/2) = x;
% 
% paddedimage = paddedimage.*prefilter_2D;
% 
% KFFTimag=fftshift(fft2(fftshift(paddedimage)));
% 
% t=1:0.0001:N;
% y=zeros(1,length(t));
% y(320000:330000)=1;
% f=fftshift(fft(y));
% f=f+1e-20+1e-20*i;
% data = giveNUFFT2D(KFFTimag,fn1,kloc,J,K,N,Ofactor);
% data1=data;
% for l=1:1:length(kloc)
%     kx=real(kloc(l))-32.5;
%     ky=imag(kloc(l))-32.5;
%     d=sqrt(kx^2+ky^2);
%     
%     data1(l)=f(floor(abs(10000*d+315000)))*data(l);
% end
%     end
% 
% 
% 
% function [idata]=INUFFT2_v(data,fn1,kloc,J,K,N,Ofactor,dcf,prefilter_2D);
% t=1:0.0001:N;
% y=zeros(1,length(t));
% y(320000:330000)=1;
% f=fftshift(fft(y));
% f=f+1e-20+1e-20*i;
% data1=data;
% for l=1:1:length(kloc)
%     kx=real(kloc(l))-32.5;
%     ky=imag(kloc(l))-32.5;
%     d=sqrt(kx^2+ky^2);
%     
%     data1(l)=data(l)*(f(floor(abs(10000*d+315000))))';
% end
% gdata = agiveNUFFT2D(data1,fn1,kloc,J,K,N,Ofactor,dcf);
% 
% 
% 
% gdata=reshape(gdata,K,K);
% 
% 
% gdata=fftshift(ifft2(fftshift(gdata)));
% 
% 
% gdata=gdata.*prefilter_2D;
% 
% 
% diff=K-N;
% idata=gdata(diff/2+1:end-diff/2,diff/2+1:end-diff/2); 
% 
% 
% end
% 
% 
% % 
% % function data=NUFFT2(x,K,diff,fn1,kloc,J,N,Ofactor,prefilter_2D)
% % paddedimage = zeros(K,K);
% % 
% % paddedimage(diff/2+1:end-diff/2,diff/2+1:end-diff/2) = x;
% % 
% % paddedimage = paddedimage.*prefilter_2D;
% % 
% % KFFTimag=fftshift(fft2(fftshift(paddedimage)));
% % 
% % 
% % data = giveNUFFT2D(KFFTimag,fn1,kloc,J,K,N,Ofactor);
% % 
% % end
% % 
% % function idata=INUFFT2(data,fn1,kloc,J,K,N,Ofactor,dcf,prefilter_2D)
% % gdata = agiveNUFFT2D(data,fn1,kloc,J,K,N,Ofactor,dcf);
% % 
% % 
% % 
% % gdata=reshape(gdata,K,K);
% % 
% % 
% gdata=fftshift(ifft2(fftshift(gdata)));
% 
% 
% gdata=gdata.*prefilter_2D;
% 
% 
% diff=K-N;
% idata=gdata(diff/2+1:end-diff/2,diff/2+1:end-diff/2); 
% end
