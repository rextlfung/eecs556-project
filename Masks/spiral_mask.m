function mask_spir = spiral_mask(nturns, rad, fov)
% nturns = number of turns (integer)
% rad = max radius of spiral
% fov = number of voxels in [nx, ny]

% given values
pos = [0 0 ;    % startpoint
       0 rad ] ;  % endpoint
%nturns = 50 ;    % number of turns (integer value)

% engine
dp = diff(pos,1,1) ;
R = hypot(dp(1), dp(2)) ;
phi0 = atan2(dp(2), dp(1)) ;

phi = linspace(0, nturns*2*pi, 200000) ; % 10000 = resolution
r = linspace(0, R, numel(phi)) ;
x = pos(1,1) + r .* cos(phi + phi0) ;
y = pos(1,2) + r  .* sin(phi + phi0) ;

%plot(x,y,'b-',pos(:,1),pos(:,2),'ro-') ; % nturns crossings, including end point

mask_spir = zeros(fov);
for l = 1:length(x)
    mask_spir(round(x(l))+fov(1)/2,round(y(l))+fov(2)/2) = 1;
end
mask_spir(128-10:128+10,128-10:128+10) = 1;

mask_spir = logical(mask_spir);