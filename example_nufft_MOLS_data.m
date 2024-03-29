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

%% Load MOLS data and get undersampled k-space from a slice
load MOLS_x_k.mat; % I, X, kloc_centered
img = permute(I, [2 1]);
img = double(img);
ksp = ifftshift(fft2(fftshift(img))); % Make "ground truth" k-space from img

% Turn raw k-space data and (kx,ky) location vectors into columns
ksp_mols = X.'; kloc_centered = kloc_centered.';

%% map non-cartesian sampled data onto high-res grid for plotting
kxx = imag(kloc_centered);
kyy = real(kloc_centered);
ksp_us = zeros(2*size(ksp));

% map k-space locations onto high-res grid indices
kxx_mapped = round(rescale(kxx, 1, size(ksp_us,1)));
kyy_mapped = round(rescale(kyy, 1, size(ksp_us,2)));

% allocate raw k-space data to closest point on grid
for n = 1:numel(ksp_mols)
    ksp_us(kxx_mapped(n), kyy_mapped(n)) = ksp_mols(n);
end

%% create NUFFT structure
N = size(ksp);
J = [2 2]; % interpolation neighborhood
K = N*2; % two-times oversampling
om = [kxx kyy];	% 'frequencies' are locations here!

%% NUFFT magic
st = nufft_init(om, N, J, K, N/2, 'minmax:kb');
weights = ksp_mols;
[pattern, Xk] = nufft_adj_modified(weights, st);
Xk = ifftshift(Xk);

%% Viz
close all;
figure; tiledlayout(2,3,'TileSpacing','tight')
nexttile; im(log(abs(ksp))); title('Ground truth k-space');
nexttile; im(log(abs(ksp_us))); title('Undersampled (zero-inserted) k-space');
nexttile; im(log(abs(Xk))); title('Interpolated k-space (upsampled)');
nexttile; im(img); title('Ground truth image')
nexttile; im(ifftshift(ifft2(fftshift(ksp_us)))); title('Zero-insertion recon');
nexttile; im(pattern); title('NUFFT recon')

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
	pattern = nufft_adj(weights, st);

	pl(2+choice)
	im(pattern, 'pattern')
end
