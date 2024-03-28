function mask_rad = radial_mask(r,n,fov)
% r = radius of spokes
% n = number of spokes
% fov = number of voxels in [nx, ny]

% r radius of spoke
% n = 50 % Number of spokes
rad_space = pi/n;
points = -r:r;
sp_x = zeros(2*r+1,n);
sp_y = zeros(2*r+1,n);
for i = 0:n-1
    sp_x(:,i+1) = points*cos(rad_space*i);
    sp_y(:,i+1) = points*sin(rad_space*i);
end

mask_rad = zeros(fov);

for j = 1:n
    skx = round(sp_x(:,j));
    sky = round(sp_y(:,j));
    for k = 1:2*r+1
    mask_rad(skx(k)+fov(1)/2,sky(k)+fov(2)/2) = 1;
    end
end

mask_rad = logical(mask_rad);