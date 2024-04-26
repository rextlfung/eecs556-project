% example_nufft1_antenna.m applied to mristack data
%
% This m-file is an example of applying the NUFFT method to compute
% the far-field beam pattern of a 2D array of nonuniformly spaced
% (point?) antenna elements.
% In the nomenclature of the 1999 SIAM J Sci. Comput. paper by Nguyen and Liu,
% this is "Problem 1."
%
% This provides an example of how to compute "the Fourier transform"
% of nonuniformly spaced data.
%
% The user should think very carefully about whether it is truly meaningful
% to compute "the FT" of unequally spaced data first.  It certainly is
% meaningful for attenna beam patterns, but may not always be what is really
% useful for all applications with nonuniformly sampled data.
%
% Copyright 2007, Jeff Fessler, University of Michigan

% Modified 2024, Rex Fung, University of Michigan
% This m-file is an example of applying the NUFFT method to compute
% the image of a 2D array of nonuniformly sampled k-space of matlab's
% built-in mristack dataset
clear; close all;
addpath('Masks/');

%% Options
use_dcf = 1; % Whether or not to do density compensation
us_factor = 16; % Upsampling factor when gridding for zero-insertion or MIRT

%% Load MOLS data and get undersampled k-space from a slice
load MOLS_x_k.mat; % I, X, kloc_centered
img = permute(I, [2 1]);
img = double(img);
ksp = ifftshift(fft2(fftshift(img))); % Make "ground truth" k-space from img

% Turn raw k-space data and (kx,ky) location vectors into columns
ksp_mols = X.'; kloc_centered = kloc_centered.';

%% map non-cartesian sampled data onto high-res grid for plotting
ksp_us = zeros(us_factor*size(ksp));

kxx = rescale(imag(kloc_centered), -pi, pi);
kyy = rescale(real(kloc_centered), -pi, pi);

% map k-space locations onto high-res grid indices
kxx_mapped = round(rescale(kxx, 1, size(ksp_us,1)));
kyy_mapped = round(rescale(kyy, 1, size(ksp_us,2)));

% allocate raw k-space data to closest point on grid
for n = 1:numel(ksp_mols)
    ksp_us(kxx_mapped(n), kyy_mapped(n)) = ksp_mols(n) * (0.3 + abs(kloc_centered(n)));
end
img_zi = ifftshift(ifft2(fftshift(ksp_us)));

% Extract central portion from upsampled (zero-inserted) image
range_x = ((-size(img,1)/2 + 1):size(img,1)/2) + size(img_zi,1)/2;
range_y = ((-size(img,2)/2 + 1):size(img,2)/2) + size(img_zi,2)/2;
img_zi = img_zi(range_x, range_y);
%% create NUFFT structure
N = size(ksp);
J = [5 5]; % interpolation neighborhood
K = N*us_factor; % overample by us_factor times
om = [kxx kyy];	% 'frequencies' are locations here!

%% NUFFT magic
st = nufft_init(om, N, J, K, N/2, 'minmax:kb');
if use_dcf
    dcf = (0.015 + sqrt(kxx.^2 + kyy.^2));
else
    dcf = ones(size(kxx)); % Make dcf useless to highlight its importance
end
weights = ksp_mols .* dcf;
% weights = ksp_mols;
[img_MIRT, Xk] = nufft_adj_modified(weights, st);
Xk = ifftshift(Xk);

%% WOLS stuff
addpath('NUFFT_MOLS/NUFFT') 
addpath('NUFFT_MOLS/DFT')
addpath('NUFFT_MOLS/Optimization')

% load the original image
I =double(rgb2gray(imread('NUFFT_MOLS/Data/knee_mri_picture_z.png')));
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

% Compute and set up NUFFT forward and backward models: one time task
% These forward models are saved to file to reduce repeated computation
%----------------------------------------------------------------------------

% NUFFT parameters
K=130;
J=6;
Ofactor=151;
Order=2;
diff=K-N;

path = ['NUFFT_MOLS/Precomputed_Kernels/',num2str(K),'_',num2str(N),'_',num2str(Ofactor),'_',num2str(Order),'_',num2str(J),'_'];

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

[img_WOLS, Xk_WOLS] = INUFFT2_Symmetric(X,fn_WOLS,kloc_centered,J,K,N,Ofactor,dcf,prefilter_WOLS); %apply A_adj_WOLS to the radial k-measurements and INUFFT an image
img_WOLS = img_WOLS/norm(img_WOLS);
img_WOLS = img_WOLS.';

%% Performance metrics
NRMSE_zi = norm(rescale(real(img_zi)) - rescale(img), 'fro')/norm(rescale(img), 'fro');
PSNR_zi = psnr(rescale(real(img_zi)), rescale(img));
SSIM_zi = ssim(rescale(real(img_zi)), rescale(img));

NRMSE_MIRT = norm(rescale(real(img_MIRT)) - rescale(img), 'fro')/norm(rescale(img), 'fro');
PSNR_MIRT = psnr(rescale(real(img_MIRT)), rescale(img));
SSIM_MIRT = ssim(rescale(real(img_MIRT)), rescale(img));

NRMSE_WOLS = norm(rescale(real(img_WOLS)) - rescale(img), 'fro')/norm(rescale(img), 'fro');
PSNR_WOLS = psnr(rescale(real(img_WOLS)), rescale(img));
SSIM_WOLS = ssim(rescale(real(img_WOLS)), rescale(img));

%% Viz
close all;
figure; tiledlayout(2,4,'TileSpacing','tight')

nexttile; im(log(abs(ksp))); title('Ground truth k-space'); colorbar
nexttile; im(log(abs(ksp_us))); title('Gridded + zero-inserted k-space (upsampled)'); colorbar
nexttile; im(log(abs(Xk))); title('MIRT interpolated k-space (upsampled)'); colorbar
nexttile; im(log(abs(Xk_WOLS))); title('WOLS interpolated k-space'); colorbar

nexttile; im(img); title('Ground truth image'); colorbar
nexttile; im(img_zi); title('Gridding + zero-insertion recon'); colorbar
xlabel(sprintf('Gridding + zero-insertion: NRMSE: %f, PSNR: %f, SSIM: %f\n', NRMSE_zi, PSNR_zi, SSIM_zi));
nexttile; im(img_MIRT); title('NUFFT recon'); colorbar
xlabel(sprintf('MIRT: NRMSE: %f, PSNR: %f, SSIM: %f\n', NRMSE_MIRT, PSNR_MIRT, SSIM_MIRT));
nexttile; im(img_WOLS); title('WOLS recon'); colorbar
xlabel(sprintf('WOLS: NRMSE: %f, PSNR: %f, SSIM: %f\n', NRMSE_WOLS, PSNR_WOLS, SSIM_WOLS));

return
%% Original code by Jeff
% antenna element locations in 2d plane
%
clf, pl = @(p) subplot(220+p);
for choice=1:2
	tmp = [0:40]'/41 * 2 * pi;
	yc = sin(choice*tmp); % funny bowtie pattern for illustration
	xc = cos(tmp); clear tmp
	pl(choice), plot(xc, yc, 'o'), axis square
	xlabel 'x', ylabel 'x', title 'antenna locations'

	% create NUFFT structure
	N = [1 1]*2^8;
	J = [5 5];	% interpolation neighborhood
	K = N*2;	% two-times oversampling
	om = [xc yc];	% 'frequencies' are locations here!

	% the following line probably needs to be changed
	% to get the proper scaling/dimensions in the pattern
	% but for now i just make it fill up [-pi/20,pi/20]
	% in hopes of getting a nice 'picture' of pattern
	om = (pi/20) * om / max(om(:));
	st = nufft_init(om, N, J, K, N/2, 'minmax:kb');

	weights = ones(size(xc)); % equal weights on each element; could change

	% call the *adjoint* NUFFT; this is what does "the FT of unequal data"
	img_MIRT = nufft_adj(weights, st);

	pl(2+choice)
	im(img_MIRT, 'pattern')
end
