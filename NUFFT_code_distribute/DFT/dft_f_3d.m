function X=dft_f_3d(p3d,kloc,N,Nz,t,kz,f)

kx=real(kloc);
ky=imag(kloc);
X=zeros(size(kx));
for x=-N/2:1:N/2-1;
    for y=-N/2:1:N/2-1;
        for z=-Nz/2:1:Nz/2-1;
            X=X+(p3d(x+N/2+1,y+N/2+1,z+Nz/2+1).*exp(-1j*2*pi*t*f(x+N/2+1,y+N/2+1,z+Nz/2+1)))*exp(-2*pi*1i*kx*x/N).*exp(-2*pi*1i*ky*y/N).*exp(-2*pi*1i*kz*z/Nz);
        end
    end
end