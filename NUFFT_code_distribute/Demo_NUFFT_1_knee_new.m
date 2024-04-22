% Demo to demonstrate the benefit of Optimal least square NUFFT schemes 
%
% [1] M. Jacob,"Optimized least square non uniform fast Fourier transform (OLS-NUFFT)" , 
%     IEEE Transactions of Signal Processing, vol. 57, issue 6, pp. 2165-2177, Feb 2009 
% [2] Z. Yang and M. Jacob,"Mean square optimal NUFFT approximation for efficient non-Cartesian 
%     MRI reconstruction, JMRI, submitted
%  
%  This demo program sets up and compares the different NUFFT schemes in
%  the context of non-Cartesian MRI. 
%-------------------------------------------+
% Author: 
% Zhili Yang @ University of Rochester, Mathews Jacob @ University of Iowa
% Email: zhyang@ece.rochester.edu
%        mathews-jacob@uiowa.edu
%-------------------------------------------+

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
figure;imagesc(abs(I));axis off;colormap(gray);title('original image')
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
%%
A_WOLS = @ (z) NUFFT2_Symmetric(z,K,diff,fn_WOLS,kloc_centered,J,N,Ofactor,prefilter_WOLS);
At_WOLS = @ (z) INUFFT2_Symmetric(z,fn_WOLS,kloc_centered,J,K,N,Ofactor,dcf,prefilter_WOLS);

%% MOLS-U
%----------------

if(~exist([path,'MOLS.mat']))
    Energy_distribution = ones(1,N)';
    [prefilter_MOLS,fn_MOLS,error]=NUFFT_Assymetric(J,N,Ofactor,K,Order,Energy_distribution);                                   
    save([path,'MOLS.mat'], 'prefilter_MOLS','fn_MOLS');
else
    load([[path,'MOLS.mat']]);
end
%%
A_MOLS = @ (z)NUFFT2D_general(z,K,fn_MOLS,fn_MOLS,kloc_centered,J,N,Ofactor,prefilter_MOLS);
At_MOLS = @ (z)INUFFT2D_general(z,fn_MOLS,fn_MOLS,kloc_centered,J,K,N,Ofactor,wi,prefilter_MOLS);

%% MOLS-G
% MOLS w/ gaussian energy distribution
%----------------
path_new = ['Precomputed_Kernels/',num2str(K),'_',num2str(N),'_',num2str(Ofactor),'_',num2str(Order),'_',num2str(J),'_'];
if(~exist([path,'MOLS_G.mat']))
    mu=N/2;
    sigma=30;
    gauss=exp(-((1:N)-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
    gauss(1)=0;
    gauss(end)=0;
    Energy_distribution = gauss';
    [prefilter_MOLS_G,fn_MOLS_G,error]=NUFFT_Assymetric(J,N,Ofactor,K,Order,Energy_distribution);                                   
    save([path,'MOLS_G.mat'], 'prefilter_MOLS_G','fn_MOLS_G');
else
    load([[path,'MOLS_G.mat']]);
end

A_MOLS_G = @ (z)NUFFT2D_general(z,K,fn_MOLS_G,fn_MOLS_G,kloc_centered,J,N,Ofactor,prefilter_MOLS_G);
At_MOLS_G = @ (z)INUFFT2D_general(z,fn_MOLS_G,fn_MOLS_G,kloc_centered,J,K,N,Ofactor,wi,prefilter_MOLS_G);

%% MOLS-T
% MOLS w/ "ground-truth" energy distribution
%----------------
path_new = ['Precomputed_Kernels/',num2str(K),'_',num2str(N),'_',num2str(Ofactor),'_',num2str(Order),'_',num2str(J),'_'];
if(~exist([path,'MOLS_T.mat']))
    FT_I=fftshift(fft(ifftshift(I(:)))); %(128^2,1)
    FT_I_ds=FT_I(length(FT_I)/2-64+1:length(FT_I)/2+64); %(128,1) lower-pass filter
    I_ds=ift(ifftshift(FT_I_ds)); %down-sample to (128,1)
    Energy_distribution = abs(I_ds);
    [prefilter_MOLS_T,fn_MOLS_T,error]=NUFFT_Assymetric(J,N,Ofactor,K,Order,Energy_distribution);                                   
    save([path,'MOLS_T.mat'], 'prefilter_MOLS_T','fn_MOLS_T');
else
    load([[path,'MOLS_T.mat']]);
end

A_MOLS_T = @ (z)NUFFT2D_general(z,K,fn_MOLS_T,fn_MOLS_T,kloc_centered,J,N,Ofactor,prefilter_MOLS_T);
At_MOLS_T = @ (z)INUFFT2D_general(z,fn_MOLS_T,fn_MOLS_T,kloc_centered,J,K,N,Ofactor,dcf,prefilter_MOLS_T);



%% LS-KB OLS NUFFT
%-----------------------

diff=K-N;
[prefilter_ls,fn_ls] = giveLSInterpolator(N,K,Ofactor,J);
%prefilter_ls_2D=prefilter_ls'*prefilter_ls;
A_LS = @ (z) NUFFT2_Symmetric(z,K,diff,fn_ls,kloc_centered,J,N,Ofactor,prefilter_ls);
At_LS = @ (z) INUFFT2_Symmetric(z,fn_ls,kloc_centered,J,K,N,Ofactor,ones(size(kloc_centered)),prefilter_ls);

% LS-KB with K=2N OLS NUFFT
%-----------------------
Kbig=256;
diffBig=Kbig-N;
Jbig=4;
[prefilter_ls_2N,fn_ls_2N] = giveLSInterpolator(N,Kbig,Ofactor,Jbig);
A_LS_2N = @ (z) NUFFT2_Symmetric(z,Kbig,diffBig,fn_ls_2N,kloc_centered,Jbig,N,Ofactor,prefilter_ls_2N);
At_LS_2N = @ (z) INUFFT2_Symmetric(z,fn_ls_2N,kloc_centered,Jbig,Kbig,N,Ofactor,ones(size(kloc_centered)),prefilter_ls_2N);



%% Solve  inverse problems

% Optimization parameters
%------------------------

opts.mu = 5e-5;            % regularization parameter
opts.beta=1;               % Initial beta value; scalar used for penalty term
opts.betarate = 1.5;         % Increment use for beta
opts.ContThreshold = 1e-2; % Threshold for updating beta; increment beta when error saturates
opts.iter = 5;            % no of iterations
opts.outiter = 10;
opts.CGiterations = 55;    % Number of maximum CG iterations 
opts.Threshold = 1e-9;     % Exit condition


%
% MOLS-U
%----------------
x_init = At_MOLS(X);
[x_MOLS,cost,datacons, tvnorm] = TVrecon(A_MOLS,At_MOLS,x_init,(X),opts);

figure(1);imshow(real(x_MOLS),[0,1],'InitialMagnification','fit');axis off
title(['MOLS-U recon image;  SNR=',num2str(SNR_2D(I,x_MOLS))]);
figure(2);imshow(abs(x_MOLS-I),[0,0.1],'InitialMagnification','fit');axis off
title(['MOLS error image']);

% WOLS-U
%----------------

x_init = At_WOLS(X);
[x_WOLS,cost,datacons, tvnorm] = TVrecon(A_WOLS,At_WOLS,x_init,(X),opts);

figure(3);imshow(abs(x_WOLS),[0,1],'InitialMagnification','fit');axis off
title(['WOLS recon image;  SNR=',num2str(SNR_2D(I,x_WOLS))]);
figure(4);imshow(abs(x_WOLS-I),[0,0.1],'InitialMagnification','fit');axis off
title(['WOLS error image']);


% Kaiser Bessel
%----------------

x_init = At_LS(X);
[x_KB,cost,datacons, tvnorm] = TVrecon(A_LS,At_LS,x_init,(X),opts);

figure(5);imshow(abs(x_KB),[0,1],'InitialMagnification','fit');axis off
title(['KB recon image;  SNR=',num2str(SNR_2D(I,x_KB))]);
figure(6);imshow(abs(x_KB-I),[0,0.1],'InitialMagnification','fit');axis off
title(['KB error image']);

%
% Kaiser Bessel: K=2N 
%----------------

x_init = At_LS_2N(X);
[x_KB_2N,cost,datacons, tvnorm] = TVrecon(A_LS_2N,At_LS_2N,x_init,(X),opts);

figure(7);imshow(abs(x_KB_2N),[0,1],'InitialMagnification','fit');axis off
title(['K=2N recon image;  SNR=',num2str(SNR_2D(I,x_KB_2N))]);
figure(8);imshow(abs(x_KB_2N-I),[0,0.1],'InitialMagnification','fit');axis off
title(['K=2N error image']);
