function [Kernel,calcweight,error] = calcKernelDiscretemod_fm(DFTMtx, fn,K,N,Ofactor, Order,H)

vector = [-K*Ofactor/2:1:K*Ofactor/2-1];
abeta = transpose(fftshift(fft(fftshift(Bspline(vector,2*(Order)-1)))));

% FT from -2*pi*Ofactor/2 to 2*pi*Ofactor/2
vector = 2*pi*[-K*Ofactor/2:1:K*Ofactor/2-1]/K;
BsplineFourier = sinc_new(vector/(2*Ofactor)).^(2*(Order))';

index = (Ofactor-1)*K/2;
bbeta = abeta;
bbeta(index+1:index+K) = abeta(index+1:index+K) - BsplineFourier(index+1:index+K);


    FN = (DFTMtx*fn);
    
    FN_den = abs(FN).^2.*abeta;
    FN_num = abs(FN).^2.*bbeta;
    
    Den=zeros(K,1);Num=zeros(K,1);
    for i=0:Ofactor-1,
        Den = Den + FN_den(i*K+1:(i+1)*K);
        Num = Num + FN_num(i*K+1:(i+1)*K);
    end
    Kernel =(Num./Den);
    Kernel([1:K/2-N/2,K/2+N/2+2:end])=0;
    
    H2=zeros(K,1);
    H2(1+K/2-N/2:end-K/2+N/2)=H;
    error = abs(mean(H2.*(Kernel)));
    calcweight = (H2./(Den))';
     
    calcweight([1:K/2-N/2,K/2+N/2+2:end])=0;
    CentralIndex = (Ofactor-1)/2*K;
    calcweight = repmat(calcweight,1,Ofactor);
    calcweight(CentralIndex+1:CentralIndex+K)=0;