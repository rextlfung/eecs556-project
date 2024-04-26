% Demo to show the worst case optimal least square (WOLS) NUFFT scheme by
% Matt. Jacobs @UofIowa
% Some mex functions of the original codes are broken so I keep the (very
% liitle) part that can still run;
% Modified 2024 Yongli He, U of M
%
%%
clear;
clc;
addpath('./NUFFT') 
addpath('./DFT')
addpath('./Optimization')

%%
% load the original image
I =double(rgb2gray(imread('Data/knee_mri_picture_z.png')));
I=I./max(I(:));
%figure;imagesc(abs(I));axis off;colormap(gray);title('original image')
N=size(I,1);

% generate a radial trajectory with 128 lines.
kloc_onesided=getpolar(128,N);
kloc_centered=kloc_onesided-N/2-N/2*1i-1-1i;

% Compute the exact Fourier samples on the radial trajectory. 
N=size(I,1);
M=length(kloc_centered);
[A1,B1]=creatA1(kloc_centered,N);
[A2,B2]=creatA(kloc_centered,N);

% Simulate the k-space measurements on the radial trajectory using DTFT. 
X=NFT_n(I,N,A1,B1,M);
% NUFFT parameters
[K Ofactor J ]=deal(N+2,101,4);
diff=K-N;
Order=2;

%% Compute and set up NUFFT forward and backward models: one time task
% These forward models are saved to file to reduce repeated computation
%----------------------------------------------------------------------------

% NUFFT parameters
K=130;
J=6;
Ofactor=151;
Order=2;
diff=K-N;
dcf=(0.3 + abs(kloc_centered));

path = ['Precomputed_Kernels/',num2str(K),'_',num2str(N),'_',num2str(Ofactor),'_',num2str(Order),'_',num2str(J),'_'];

% Worst case OLS
%----------------

if(~exist([path,'WOLS.mat']))
    [prefilter_WOLS,fn_WOLS]=NUFFT_symmetric(J,N,Ofactor,K,Order);
    save([path,'WOLS.mat'], 'prefilter_WOLS','fn_WOLS');
else
    load([[path,'WOLS.mat']]);
end

A_WOLS = @ (z) NUFFT2_Symmetric(z,K,diff,fn_WOLS,kloc_centered,J,N,Ofactor,prefilter_WOLS);
At_WOLS = @ (z) INUFFT2_Symmetric(z,fn_WOLS,kloc_centered,J,K,N,Ofactor,dcf,prefilter_WOLS);

x_WOLS = At_WOLS(X); %apply A_adj_WOLS to the radial k-measurements and INUFFT an image

%% MIRT NUFFT
% Turn raw k-space data and (kx,ky) location vectors into columns
ksp_mols = X.'; kloc_centered = kloc_centered.';
kxx = rescale(real(kloc_centered), -pi, pi);
kyy = rescale(imag(kloc_centered), -pi, pi);
% create NUFFT structure
N = [size(I)];
J = [5 5]; % interpolation neighborhood
K = N*2; % two-times oversampling
om = [kxx, kyy]; % 'frequencies' are locations here!

% Desnity compensated weights
weights = ksp_mols .* dcf.';
% weights = ksp_mols;

% NUFFT magic
st = nufft_init(om, N, J, K, N/2, 'minmax:kb');
x_MIRT = nufft_adj(weights,st);

% COMPUTE METRICS
PSNR_MIRT = psnr(rescale(real(x_MIRT)), rescale(I));
SSIM_MIRT = ssim(rescale(real(x_MIRT)), rescale(I));

PSNR_WOLS = psnr(rescale(real(x_WOLS)), rescale(I));
SSIM_WOLS = ssim(rescale(real(x_WOLS)), rescale(I));

% DISPLAY
figure; tiledlayout(1,3,'TileSpacing','tight')
nexttile;im(abs(I')); title('Original image'); colorbar;
nexttile;im(abs(x_MIRT')); title('NUFFT_MIRT'); colorbar;
xlabel(sprintf("PSNR: %f, SSIM: %f", PSNR_MIRT, SSIM_MIRT))
nexttile;im(abs(x_WOLS')); title('NUFFT_WOLS'); colorbar;
xlabel(sprintf("PSNR: %f, SSIM: %f", PSNR_WOLS, SSIM_WOLS))