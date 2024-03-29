% Demo to show the worst case optimal least square (WOLS) NUFFT scheme by
% Matt. Jacobs @UofIowa
% Some mex functions of the original codes are broken so I keep the (very
% liitle) part that can still run;
% Modified 2024 Yongli He, U of M
%
%%
clear all;
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
dcf=ones(1,length(kloc_centered));

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
At_WOLS = @ (z) INUFFT2_Symmetric(z,fn_WOLS,kloc_centered,J,K,N,Ofactor,ones(size(kloc_centered)),prefilter_WOLS);

x_WOLS = At_WOLS(X); %apply A_adj_WOLS to the radial k-measurements and INUFFT an image

figure; tiledlayout(1,2,'TileSpacing','tight')
nexttile;im(abs(I')); title('Original image')
nexttile;im(abs(x_WOLS')); title('NUFFT_WOLS')
