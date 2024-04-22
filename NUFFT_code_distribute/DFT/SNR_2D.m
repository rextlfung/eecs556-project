%Signal to noise ratio

function SNR=SNR_2D(u0,u)
%  u0: orginal image
%  u: noised image
%  (nb,na): size 
if size(u0)~=size(u)
    disp('They are not the same size')
    exit
end

[nb,na]=size(u0);
MSE=sum((abs(u(:))-abs(u0(:))).^2);
SNR=10.*log10(sum(abs(u0(:)).^2)/MSE);

return
