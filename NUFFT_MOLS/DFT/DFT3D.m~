function X=DFT3D(p3d,kloc,N)
kz=squeeze(kloc(2,:));
kx=real(squeeze(kloc(1,:)));
ky=imag(squeeze(kloc(1,:)));
X=zeros(size(kx));
for x=-N/2:1:N/2-1;
    for y=-N/2:1:N/2-1;
        for z=-N/2:1:N/2-1;
            X=X+p3d(x+N/2+1,y+N/2+1,z+N/2+1)*exp(-2*pi*i*kx*x/N).*exp(-2*pi*i*ky*y/N).*exp(-2*pi*i*kz*z/N);
        end
    end
end