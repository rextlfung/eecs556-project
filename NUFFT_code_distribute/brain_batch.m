% read 20 T1w brain mri structural images in nifti format.
% find the prior distribution of this dataset by averaging the
% lowpass-filtered spectra

data_root='/Users/yonglismac/Documents/WIN24/556Project/brainmri_nifti/train';

files = dir(fullfile(data_root, '*')); %first two elements '.' and '..' should be excluded
I_batch=zeros(256,256,numel(files)-3);
%%
for i=1:numel(files)-3
    i
    imFile=fullfile(data_root,files(3+i).name);
    metadata=niftiinfo(imFile);
    vol=niftiread(metadata);
    [nx,ny,nz]=size(vol);
    I_batch(1:nx,1:ny,i)=vol(:,:,round(nz/2));
end

idx=[1,2,3,4,5,7,8,11,12,13,14,15,16,17,19];
I_batch_15=I_batch(:,:,idx);
% I_batch_15=I_batch;

batchsize=size(I_batch_15,3);
% image preprocessing
for i=1:size(I_batch_15,3)
    I_batch_15(:,:,i)=I_batch_15(:,:,i)-min(I_batch_15(:,:,i),[],'all');% shift to [0,X]
    I_batch_15(:,:,i)=I_batch_15(:,:,i)./max(I_batch_15(:,:,i),[],'all'); %normalize to [0,1]
end

%% Generate k samples
% generate a radial trajectory with 128 lines.
N=size(I_batch_15,1);
kloc_onesided=getpolar(128,N);
kloc_centered=kloc_onesided-N/2-N/2*1i-1-1i;

% Compute the exact Fourier samples on the radial trajectory. 
M=length(kloc_centered);
[A1,B1]=creatA1(kloc_centered,N);
[A2,B2]=creatA(kloc_centered,N);

% Simulate the k-space measurements on the radial trajectory using DTFT. 
X_batch=zeros(M,batchsize);
for i=1:batchsize
    i
    X_batch(:,i)=NFT_n(I_batch_15(:,:,i),N,A1,B1,M);
end
%% NUFFT parameters
K=N+4;
J=6;
Ofactor=151;
Order=2;
diff=K-N;
wi = (0.3 + abs(kloc_centered));

path = ['Precomputed_Kernels/',num2str(K),'_',num2str(N),'_',num2str(Ofactor),'_',num2str(Order),'_',num2str(J),'_'];

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
    sigma=30*round(N/128);
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
im(abs(reshape(I_batch_15(:,:,1),[],1)));title('I1');xlabel('n');ylabel('a.u.'); %one example image signal
figure;im(abs(reshape(I_batch_15(:,:,3),[],1)));title('I2');xlabel('n');ylabel('a.u.') % another example signal
I_avg=mean(I_batch_15,3); % average image
figure;im(I_avg); %display average image 
figure;im(abs(I_avg(:)));title('Average Signal');xlabel('n');ylabel('a.u.') %display average signal

path_new = ['Precomputed_Kernels/',num2str(K),'_',num2str(N),'_',num2str(Ofactor),'_',num2str(Order),'_',num2str(J),'_'];
if(~exist([path,'MOLS_T.mat']))
    FT_I=fftshift(fft(ifftshift(I_avg(:)))); %(N^2,1)
    FT_I_ds=FT_I(length(FT_I)/2-N/2+1:length(FT_I)/2+N/2); %(N,1) lower-pass filter
    I_ds=ift(ifftshift(FT_I_ds)); %down-sample to (N,1)
    Energy_distribution = abs(I_ds);
    [prefilter_MOLS_T,fn_MOLS_T,error]=NUFFT_Assymetric(J,N,Ofactor,K,Order,Energy_distribution);                                   
    save([path,'MOLS_T.mat'], 'prefilter_MOLS_T','fn_MOLS_T');
else
    load([[path,'MOLS_T.mat']]);
end

A_MOLS_T = @ (z)NUFFT2D_general(z,K,fn_MOLS_T,fn_MOLS_T,kloc_centered,J,N,Ofactor,prefilter_MOLS_T);
At_MOLS_T = @ (z)INUFFT2D_general(z,fn_MOLS_T,fn_MOLS_T,kloc_centered,J,K,N,Ofactor,wi,prefilter_MOLS_T);

%% Quantative metrics
xinit_U=zeros(size(I_batch_15));
xinit_G=zeros(size(I_batch_15));
xinit_T=zeros(size(I_batch_15));
NRMSE_U=zeros(1,batchsize); NRMSE_G=zeros(1,batchsize); NRMSE_T=zeros(1,batchsize);
PSNR_U=zeros(1,batchsize); PSNR_G=zeros(1,batchsize); PSNR_T=zeros(1,batchsize); 
for i=1:batchsize
    i
    xinit_U(:,:,i)=At_MOLS(X_batch(:,i));
    xinit_G(:,:,i)=At_MOLS_G(X_batch(:,i));
    xinit_T(:,:,i)=At_MOLS_T(X_batch(:,i));

    NRMSE_U(i) = norm(rescale(real(xinit_U(:,:,i))) - rescale(I_batch_15(:,:,i)),'fro')/norm(rescale(I_batch_15(:,:,i)),'fro');
    NRMSE_G(i) = norm(rescale(real(xinit_G(:,:,i))) - rescale(I_batch_15(:,:,i)),'fro')/norm(rescale(I_batch_15(:,:,i)),'fro');
    NRMSE_T(i) = norm(rescale(real(xinit_T(:,:,i))) - rescale(I_batch_15(:,:,i)),'fro')/norm(rescale(I_batch_15(:,:,i)),'fro');

    PSNR_U(i) = psnr(rescale(real(xinit_U(:,:,i))), rescale(I_batch_15(:,:,i)));
    PSNR_G(i) = psnr(rescale(real(xinit_G(:,:,i))), rescale(I_batch_15(:,:,i)));
    PSNR_T(i) = psnr(rescale(real(xinit_T(:,:,i))), rescale(I_batch_15(:,:,i)));

end

NRMSE_U_avg=mean(NRMSE_U);
NRMSE_G_avg=mean(NRMSE_G);
NRMSE_T_avg=mean(NRMSE_T);

PSNR_U_avg=mean(PSNR_U);
PSNR_G_avg=mean(PSNR_G);
PSNR_T_avg=mean(PSNR_T);

%% testing
%% testing dataset
data_root_tst='/Users/yonglismac/Documents/WIN24/556Project/brainmri_nifti/test';

files_tst = dir(fullfile(data_root_tst, '*')); %first two elements '.' and '..' should be excluded
I_batch_tst=zeros(256,256,numel(files_tst)-2); %exclude first two './' and '../'

for i=1:numel(files_tst)-2
    i
    imFile=fullfile(data_root_tst,files_tst(2+i).name);
    metadata=niftiinfo(imFile);
    vol=niftiread(metadata);
    [nx,ny,nz]=size(vol);
    I_batch_tst(1:nx,1:ny,i)=vol(:,:,round(nz/2));
end

I_batch_tst=I_batch_tst(:,:,[1,2,4,5,6]); %take 5 images
bs_tst=size(I_batch_tst,3);
% image preprocessing
for i=1:bs_tst
    I_batch_tst(:,:,i)=I_batch_tst(:,:,i)-min(I_batch_tst(:,:,i),[],'all');% shift to [0,X]
    I_batch_tst(:,:,i)=I_batch_tst(:,:,i)./max(I_batch_tst(:,:,i),[],'all'); %normalize to [0,1]
end
%% testing images k-space
% Simulate the k-space measurements on the radial trajectory using DTFT. 

X_batch_tst=zeros(M,bs_tst);
for i=1:bs_tst
    i
    X_batch_tst(:,i)=NFT_n(I_batch_tst(:,:,i),N,A1,B1,M);
end
%% testing metrics
xinit_U_tst=zeros(size(I_batch_tst));
xinit_G_tst=zeros(size(I_batch_tst));
xinit_T_tst=zeros(size(I_batch_tst));
NRMSE_tst_U=zeros(1,bs_tst); NRMSE_tst_G=zeros(1,bs_tst); NRMSE_tst_T=zeros(1,bs_tst);
PSNR_tst_U=zeros(1,bs_tst); PSNR_tst_G=zeros(1,bs_tst); PSNR_tst_T=zeros(1,bs_tst); 
for i=1:bs_tst
    i
    xinit_U_tst(:,:,i)=At_MOLS(X_batch_tst(:,i));
    xinit_G_tst(:,:,i)=At_MOLS_G(X_batch_tst(:,i));
    xinit_T_tst(:,:,i)=At_MOLS_T(X_batch_tst(:,i));

    NRMSE_tst_U(i) = norm(rescale(real(xinit_U_tst(:,:,i))) - rescale(I_batch_tst(:,:,i)),'fro')/norm(rescale(I_batch_tst(:,:,i)),'fro');
    NRMSE_tst_G(i) = norm(rescale(real(xinit_G_tst(:,:,i))) - rescale(I_batch_tst(:,:,i)),'fro')/norm(rescale(I_batch_tst(:,:,i)),'fro');
    NRMSE_tst_T(i) = norm(rescale(real(xinit_T_tst(:,:,i))) - rescale(I_batch_tst(:,:,i)),'fro')/norm(rescale(I_batch_tst(:,:,i)),'fro');

    PSNR_tst_U(i) = psnr(rescale(real(xinit_U_tst(:,:,i))), rescale(I_batch_tst(:,:,i)));
    PSNR_tst_G(i) = psnr(rescale(real(xinit_G_tst(:,:,i))), rescale(I_batch_tst(:,:,i)));
    PSNR_tst_T(i) = psnr(rescale(real(xinit_T_tst(:,:,i))), rescale(I_batch_tst(:,:,i)));

end

NRMSE_tst_U_avg=mean(NRMSE_tst_U);
NRMSE_tst_G_avg=mean(NRMSE_tst_G);
NRMSE_tst_T_avg=mean(NRMSE_tst_T);

PSNR_tst_U_avg=mean(PSNR_tst_U);
PSNR_tst_G_avg=mean(PSNR_tst_G);
PSNR_tst_T_avg=mean(PSNR_tst_T);

