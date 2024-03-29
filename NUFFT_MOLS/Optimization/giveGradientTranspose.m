
%-------------------------------------------
% DiffOpp = giveDuchonSecondDiffOperator(X)
% Computes the Duchon differential operators
% Assumes periodic boundary conditions
%-------------------------------------------

function X = giveGradientTranspose(DiffOpp)

XX = DiffOpp{1};
YY = DiffOpp{2};


%-------------------------------------
% XX = X(n+1,:) - 2*X(n,:) + X(n-1,:)
%-------------------------------------
X = (circshift(XX,[-1,0]) - XX);

%-------------------------------------
% YY = X(:,n+1) - 2*X(:,n) + X(:,n-1)
%-------------------------------------
temp = (circshift(YY,[0,-1]) - YY);

X = X+ temp;
