function [X,Y,samp] = sample_radial(img,r,n)
% Inputs
% img = 2D image/kspace
% r = radius of spokes
% n = number of spokes

% Outputs
% X = x coordinate of sample
% Y = y coordinate of sample
% samp = interpolated values of points
fov = size(img);
rad_space = pi/n;
points = -r:r;
sp_x = zeros(2*r+1,n);
sp_y = zeros(2*r+1,n);
for i = 0:n-1
    sp_x(:,i+1) = points*cos(rad_space*i);
    sp_y(:,i+1) = points*sin(rad_space*i);
end

X = sp_x(:);
Y = sp_y(:);
[xb,yb] = meshgrid(-fov(1)/2:fov(1)/2-1);

samp = interp2(xb,yb,img,X,Y);