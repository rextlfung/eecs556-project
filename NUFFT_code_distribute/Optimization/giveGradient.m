
%-------------------------------------------
% DiffOpp = giveGradient(X)
% give the gradients at f(k+1/2,l) and f(k,l+1/2)
% Derivatives of Bspline at f(k+1/2) are (-1,1)
% Assumes periodic conditions
%-------------------------------------------

function DiffOpp = giveGradient(X)

DiffOpp = cell(2,1);



%-------------------------------------
% XX = X(n+1,:) - 2*X(n,:) + X(n-1,:)
%-------------------------------------
XX =  circshift(X,[1,0]) - X ;
%-------------------------------------
% YY = X(:,n+1) - 2*X(:,n) + X(:,n-1)
%-------------------------------------
YY = circshift(X,[0,1]) - X ;


DiffOpp{1} = XX;
DiffOpp{2} = YY;