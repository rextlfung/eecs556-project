function [A1,B1]=creatA1(kloc,N)
kx=[];
ky=[];
kx=real(kloc);
ky=imag(kloc);

n=-N/2:1:N/2-1;
m=-N/2:1:N/2-1;
%n=1:1:N;

%m=1:1:N;

A1=exp(-1j*2*pi/N*((kx')*n));
B1=exp(-1j*2*pi/N*((ky')*m));
end
