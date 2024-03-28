function [X,Y,samp] = sample_spiral(img,r,n)

% Inputs
% img = 2D image/kspace
% r = radius of outermost spiral
% n = number of turns

% Outputs
% X = x coordinate of sample
% Y = y coordinate of sample
% samp = interpolated values of points

pos = [0 0 ;    % startpoint
       0 r ] ;  % endpoint

fov = size(img);
% engine
dp = diff(pos,1,1) ;
R = hypot(dp(1), dp(2)) ;
phi0 = atan2(dp(2), dp(1)) ;

phi = linspace(0, nturns*2*pi, 200000) ; % 10000 = resolution
r = linspace(0, R, numel(phi)) ;
X = pos(1,1) + r .* cos(phi + phi0) ;
Y = pos(1,2) + r  .* sin(phi + phi0) ;

[xb,yb] = meshgrid(-fov(1)/2:fov(1)/2-1);

samp = interp2(xb,yb,img,X,Y);
