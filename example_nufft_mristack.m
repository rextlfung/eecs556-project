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

%% Load mristack and get undersampled k-space from a slice
load mristack;
mristack = permute(mristack, [2 1 3]);
img = double(mristack(:,:,8));
ksp = ifftshift(fft2(fftshift(img)));

% Test with dummy uniform sampling mask
mask = false(size(ksp));
mask(:,1:2:end) = true;

% Undersample
ksp_us = ksp .* mask;

%% Get k-space locations from data
kx = linspace(-pi, pi, 256);
ky = linspace(-pi, pi, 256);
[kxx, kyy] = ndgrid(kx, ky);

% Undersampled version
kxx_us = kxx(mask);
kyy_us = kyy(mask);

% create NUFFT structure
N = size(ksp);
J = [5 5];	% interpolation neighborhood
K = N*2;	% two-times oversampling
om = [kxx_us(:) kyy_us(:)];	% 'frequencies' are locations here!

%% NUFFT magic
st = nufft_init(om, N, J, K, N/2, 'minmax:kb');
weights = ksp(mask);
[pattern, Xk] = nufft_adj(weights(:), st);
Xk = ifftshift(Xk);

%% Viz
figure; tiledlayout('flow','TileSpacing','tight')
nexttile; im(log(abs(ksp_us))); title('Undersampled k-space');
nexttile; im(ifftshift(ifft2(fftshift(ksp_us)))); title('Zero-insertion recon');
nexttile; im(log(abs(Xk))); title('Interpolated k-space (upsampled)');
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
