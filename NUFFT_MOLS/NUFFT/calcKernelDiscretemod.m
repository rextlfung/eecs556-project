function [Kernel,calcweight,error,bbeta] = calcKernelDiscretemod(DFTMtx, fn,K,N,Ofactor, Order,scalefactors,J)

if(nargin > 6)
    Kvector = [-K/2:K/2-1];
    mask = ((Kvector >=-N/2) &(Kvector <= N/2-1));
    index = find(mask);
    tempfn = [flipud(fn(2:end));fn];
    sv = 0*Kvector;sv(index) = (scalefactors)./Ofactor;
end

vector = [-K*Ofactor/2:1:K*Ofactor/2-1];
abeta = real(fftshift(fft(fftshift(Bspline(vector,2*(Order)-1)))))';

% FT from -2*pi*Ofactor/2 to 2*pi*Ofactor/2
vector = 2*pi*[-K*Ofactor/2:1:K*Ofactor/2-1]/K;
BsplineFourier = sinc(vector/(2*Ofactor)).^(2*(Order))';

index = (Ofactor-1)*K/2;bbeta = abeta;
bbeta(index+1:index+K) = abeta(index+1:index+K) - BsplineFourier(index+1:index+K);


FN = (DFTMtx*fn);
    
    
    FN_den = abs(FN).^2.*abeta;
    FN_num = abs(FN).^2.*bbeta;
    
    Den=zeros(K,1);Num=zeros(K,1);
    for i=0:Ofactor-1,
        Den = Den + FN_den(i*K+1:(i+1)*K);
        Num = Num + FN_num(i*K+1:(i+1)*K);
    end
    Kernel = (Num./Den);
    if(nargin>6)
        centre = length(FN)/2;
        FN1 = FN.*BsplineFourier;
        svopt =(FN1(centre-K/2+1:centre+K/2)./Den)';
        Kernel = Kernel + Den.*(abs((sv-svopt)').^2);
    end
    Kernel([1:K/2-N/2,K/2+N/2+1:end])=0;
    error = sqrt(mean(Kernel.^2));
    
     calcweight = (Num./(Den.^2))';calcweight([1:K/2-N/2,K/2+N/2+2:end])=0;
     CentralIndex = (Ofactor-1)/2*K;
     calcweight = repmat(calcweight,1,Ofactor); calcweight(CentralIndex+1:CentralIndex+K)=0;