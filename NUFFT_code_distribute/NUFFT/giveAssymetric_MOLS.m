function [fn,Kernel,e]=  giveAssymetric_MOLS(J,K,N,Ofactor,Order,H)

if(mod(Ofactor,2)==0)
    display('Ofactor cannot be even');
    return;
end
if(mod(K,2)~=0 | mod(N,2)~=0)
    display('K and N should be even');
    return;
end



x = [0:1/Ofactor:ceil(J/2)];
x = x(find(x<=J/2));
centreindex = length(x);
x = [-fliplr(x(2:end)),x];

OrigWeight = ones(K,1);


Nsamples = Ofactor*K;
k = [-Nsamples/2:Nsamples/2-1]';
DFTMtx = exp(-1j*2*pi*k*x/K);

vector = [0:1:K*Ofactor/2-1,-K*Ofactor/2:1:-1];
abeta = K*Ofactor*fftshift(ifft(fftshift(Bspline(vector,3))))';
len = length(abeta);
Weight = sparse([1:len],[1:len],abeta);
B = DFTMtx'*Weight*DFTMtx;
fn = Bspline(x(centreindex:end),J-1);
fn2=fliplr(fn(2:end));
fn1= [fn2,fn];
fn=transpose(fn1);
[Kernel,current_weight,error] = calcKernelDiscretemod_fm(DFTMtx,fn,K,N,Ofactor,Order,H);
oldfn = fn;
olderror = error;
 e=error;
CentralIndex = (Ofactor-1)/2*K;
%   if(displayflag)
% %    figure(figindex)
%     subplot(2,1,1);plot(real(fn));title(strcat('Interpolator: Error = ',num2str(error)));
%     subplot(2,1,2);plot(real(Kernel));title('Kernel');
%     pause(1);
% end
%   
%%
% Start of iteration    
% -------------------    

for iter=1:100,
    weight = current_weight;
    
    len = length(weight);
    Weight = sparse([1:len],[1:len],weight);
    A = DFTMtx'*Weight*DFTMtx;
    [a,b] = eig((A),(B));
    [junk,index] = min(diag(b));

    newfn = a(:,index);
     newfn = newfn(:)./(sum((newfn)));
    if(length(oldfn)==0)
        fn = newfn;
        oldfn = fn;
        [Kernel,current_weight,error] = calcKernelDiscretemod_fm(DFTMtx,fn,K,N,Ofactor,Order,H);
        step=1;
    else
        [Kernel,current_weight,error,fn, step] = giveOptStepDiscrete_fm(DFTMtx, fn,K,N,Ofactor, olderror, oldfn, newfn,Order,H);
        oldfn = fn;
    end
    e=error;
    olderror = error;
    oldfn = fn;
  
    
    if(step==0 )
            break;
    end

end


