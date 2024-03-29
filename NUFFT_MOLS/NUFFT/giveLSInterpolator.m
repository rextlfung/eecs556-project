% function to compute LS-KB scalefactor(prefilter) and interpolator
% interpolator---LS-KB interpolators
% scalefactor---scale factors
% J---interpolator size
% N---size of image
% K---oversampled size of image
% Ofactor--oversampling factor of interpolator
% H---energy distribution of the image
function [prefilter_2D,interpolator] = giveLSInterpolator(N,K,Ofactor,J)
m=J/2;
Samples = 0:1/Ofactor:1;
k=-N/2:1:N/2-1;
alpha=K/N;
a1=pi*(2-1/alpha);
aaa=m*sqrt(a1^2-(2*pi.*(k+0.5)/K).^2);
pre=1./(K*besseli(0,(aaa)));
l = [-m:m];
D=diag(pre);
q = [];
interpolator = zeros(1,ceil(J/2)*Ofactor);
for j=1:length(Samples),
    l= -m:m;
    q = Samples(j)+l;
    mask = find(abs(q)<m);
    l = l(mask);
    q = q(mask);
   
    T = exp(-2*pi*1i*(k+0.5)'*l/K); 
    TDT = T'*D*D*T;
    TDTi = inv(TDT);

    Tt=transpose(T);
    E =exp(-2*pi*1i*(k)*Samples(j)/K);

    bb=(TDTi*Tt*D*transpose(E));
    interpolator(round((q+m)*Ofactor+1)) = bb;
end
interpolator=real(interpolator(1:end));
if mod(m*2,2)==0
    interpolator=interpolator(2:end);
end
q = [-K/2:K/2-1];
minindex = find(q==-N/2);
maxindex = find(q==N/2-1);

scalefactor = zeros(1,K);
scalefactor(minindex:maxindex) =pre;

%scalefactor = scalefactor.*interpolator(ceil(m*Ofactor));
%interpolator = interpolator./interpolator(ceil(m*Ofactor));
%plot(real(interpolator));

%scalefactor = scalefactor.*sum(interpolator(:));
interpolator = interpolator(:)./sum(interpolator);
midindex = ceil(length(interpolator)/2);
interpolator = [interpolator([midindex:end]);0];


prefilter_2D=scalefactor'*scalefactor;
kloc = 0+0*i;
scale = NUFFT2_Symmetric(ones(N),K,K-N,interpolator,kloc,J,N,Ofactor,prefilter_2D);
prefilter_2D = prefilter_2D./scale*N*N;

