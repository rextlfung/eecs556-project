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

%% Load mristack and get undersampled k-space from a slice
load mristack;
mristack = permute(mristack, [2 1 3]);
img = double(mristack(:,:,8));
img = img/norm(img);
ksp = ifftshift(fft2(fftshift(img)));

% Create sampling mask
type = 2;

switch type
    case 0 % uniform, 2x undersampled in PE direciton
        mask = false(size(ksp));
        mask(:,1:2:end) = true;
    case 1 % radial
        r = 127;
        n = 64;
        fov = size(ksp);
        mask = radial_mask(r, n, fov);
    case 2 % spiral, uniform density
        nturns = 64;
        rad = 127;
        fov = size(ksp);
        mask = spiral_mask(nturns,rad,fov);
end

% Undersample
ksp_us = ksp .* mask;
img_zi = ifftshift(ifft2(fftshift(ksp_us)));
img_zi = img_zi/norm(img_zi);

%% Get k-space locations from data
kx = linspace(-pi, pi, size(ksp,1));
ky = linspace(-pi, pi, size(ksp,2));
[kxx, kyy] = ndgrid(kx, ky);

% Undersampled version
kxx_us = kxx(mask);
kyy_us = kyy(mask);

% Density compensated weights
dcf = (pi + sqrt(kxx.^2 + kyy.^2));
ksp_weighted = ksp .* dcf;
weights = ksp_weighted(mask);

% create NUFFT structure
N = size(ksp);
J = [5 5];	% interpolation neighborhood
K = N*2;	% two-times oversampling
om = [kxx_us(:) kyy_us(:)];	% 'frequencies' are locations here!

%% NUFFT magic
st = nufft_init(om, N, J, K, N/2, 'minmax:kb');
[img_nufft, Xk] = nufft_adj_modified(weights(:), st);
img_nufft = img_nufft/norm(img_nufft);
Xk = ifftshift(Xk);

%% Performance metrics
NRMSE_zi = norm(rescale(real(img_zi)) - rescale(img),'fro')/norm(rescale(img),'fro');
PSNR_zi = psnr(rescale(real(img_zi)), rescale(img));
SSIM_zi = ssim(rescale(real(img_zi)), rescale(img));

NRMSE_MIRT = norm(rescale(real(img_nufft)) - rescale(img),'fro')/norm(rescale(img),'fro');
PSNR_MIRT = psnr(rescale(real(img_nufft)), rescale(img));
SSIM_MIRT = ssim(rescale(real(img_nufft)), rescale(img));

%% Viz
close all;
figure; tiledlayout(2,3,'TileSpacing','tight')
nexttile; im(log(abs(ksp))); title('Ground truth k-space'); colorbar
nexttile; im(log(abs(ksp_us))); title('Undersampled (zero-inserted) k-space'); colorbar
nexttile; im(log(abs(Xk))); title('Interpolated k-space (upsampled)'); colorbar
nexttile; im(img); title('Ground truth image'); colorbar
nexttile; im(img_zi); title('Zero-insertion recon'); colorbar
xlabel(sprintf('Zero-insertion: NRMSE: %f, PSNR: %f, SSIM: %f\n', NRMSE_zi, PSNR_zi, SSIM_zi));
nexttile; im(img_nufft); title('MIRT NUFFT recon'); colorbar
xlabel(sprintf('MIRT NUFFT: NRMSE: %f, PSNR: %f, SSIM: %f\n', NRMSE_MIRT, PSNR_MIRT, SSIM_MIRT));

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
	img_nufft = nufft_adj(weights, st);

	pl(2+choice)
	im(img_nufft, 'pattern')
end
