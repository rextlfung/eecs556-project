function [fn,error,Kernel]=  giveSymmetricNUFFTfn(J,K,N,Ofactor,Order)

 

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
DFTMtx = exp(-j*2*pi*k*x/K);

DFTMtx(:,centreindex+1:end) = DFTMtx(:,centreindex+1:end) + fliplr(DFTMtx(:,1:centreindex-1));
DFTMtx = DFTMtx(:,centreindex:end);

vector = [0:1:K*Ofactor/2-1,-K*Ofactor/2:1:-1];
abeta = K*Ofactor*fftshift(ifft(Bspline(vector,3)))';
len = length(abeta);
Weight = sparse([1:len],[1:len],abeta);
B = DFTMtx'*Weight*DFTMtx;

if(J<10)
    fn = Bspline(x(centreindex:end),J-1)';
else
    fn = Bspline(x(centreindex:end),9)';
end

[Kernel,current_weight,error] = calcKernelDiscretemod(DFTMtx,fn,K,N,Ofactor,Order);
oldfn = fn;
olderror = error;
 

CentralIndex = (Ofactor-1)/2*K;  
%%
% Start of iteration    
% -------------------    

for iter=1:100,
    
    
    
    weight = current_weight;
    
    len = length(weight);
    Weight = sparse([1:len],[1:len],weight);
    A = DFTMtx'*Weight*DFTMtx;
    

    [a,b] = eig(real(A),real(B));
    [junk,index] = min(diag(b));


    newfn = a(:,index);
    newfn = newfn(:)./(sum(newfn(2:end))*2+newfn(1));
    
    if(length(oldfn)==0)
        fn = newfn;
        oldfn = fn;
        [Kernel,current_weight,error] = calcKernelDiscretemod(DFTMtx,fn,K,N,Ofactor,Order);
        step=1;
    else
        [Kernel,current_weight,error,fn, step] = giveOptStepDiscrete(DFTMtx, fn,K,N,Ofactor, olderror, oldfn, newfn,Order);
        oldfn = fn;
    end


    
    olderror = error;
    oldfn = fn;
    
  
    if(step==0 )
            break;
    end

end


