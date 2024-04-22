function prefilter = givePrefilter(fn,J,K,Ofactor,Order)

% x values over the range of the interpolation coefficients
x = [0:1/Ofactor:ceil(J/2)];x = x(find(x<=J/2));
x = [-fliplr(x(2:end)),x];
Nsamples = Ofactor*K;
% K samples: ranges from -2*pi*(Ofactor-1)/2 to 2*pi*(Ofactor-1)/2
k = [-Nsamples/2:Nsamples/2-1]';
DFTMtx = exp(-1j*2*pi*k*x/K);

% DFT of the interpolation coefficients
FN = DFTMtx*fn;

% Continuous Fourier transform of the Bspline function
BsplineFourier = sinc_new(pi*k'/Nsamples).^(Order)'/Ofactor;

k = [-49*Nsamples/2:49*Nsamples/2-1]';
junk = sinc_new(pi*k'/Nsamples).^(Order)'/Ofactor;
BSPCorr = 0*BsplineFourier;
for i=0:48,
    BSPCorr = BSPCorr + abs(junk(i*Nsamples+1:i*Nsamples+Nsamples)).^2;
end

% Evaluating the denomiator
BSPCorr = BSPCorr.*abs(FN).^2;
%BSPCorr = abs(FN).^2;

Den = zeros(K,1);
for i=0:Ofactor-1,
    Den = Den + abs(BSPCorr(i*K+1:i*K+K));
end
Den = Den;

midindex = (Ofactor-1)/2;
prefilter = real(FN(midindex*K+1:midindex*K+K)).*BsplineFourier(midindex*K+1:midindex*K+K);
prefilter = prefilter./Den;

