function [A,B]=creatA(kloc,N)
kx=[];
ky=[];
kx=real(kloc);
ky=imag(kloc);
n=-N/2:1:N/2-1;
%n=1:1:N;
A=exp(j*2*pi/N*(kx')*n)*1/N;
B=exp(j*2*pi/N*(ky')*n)*1/N;
end

