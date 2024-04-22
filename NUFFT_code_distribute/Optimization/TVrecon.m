function [U,cost,datacons, tvnorm] = TVrecon(A,At,x_init,b,optsin)

opts = optsin;

% Define the forward and backward gradient operators
%----------------------------------------------------
D = @ (z) giveGradient(z);
Dt = @ (z) giveGradientTranspose(z);
U=x_init; [m n] = size(U); 

% Lam are the Lagrange multipliers. 
% refer SG Lingala et al, ISBI 2011 for details

Dfwd=D(U);
Wt{1} = 1e-18*Dfwd{1};
Wt{2} = 1e-18*Dfwd{2};
Lam = Wt;

cost=[];tvnorm=[];datacons=[];

% ================================
%  Begin Alternating Minimization
% ================================
        
for out = 1:opts.outiter,  
    sprintf('outer loop:%d',out)
    for in = double(1:opts.iter)
        sprintf('inner loop:%d',in)
        e = A(U) - b;     
        Dfwd = D(U);
        V1 = sqrt(abs(Dfwd{1}.^2)+abs(Dfwd{2}.^2));
        V1 = sum(abs(V1(:)));
      
        cost = [cost, sum(abs(e(:)).^2)   + V1*opts.mu ];
        tvnorm =[tvnorm, V1*opts.mu];
        datacons =[datacons, sum(abs(e(:)).^2)];

        Z1 = Dfwd{1} + Lam{1}/opts.beta;
        Z2 = Dfwd{2} + Lam{2}/opts.beta;
        
        % Shrinkage: update Wt
        %----------------------
        V = abs(Z1).^2+abs(Z2).^2; 
        V = sqrt(V);
        V(V==0) = 1;
        thresh_TV = (1/opts.beta);
    
        V = max(V - thresh_TV, 0)./V;
        Wt{1} = Z1.*V;
        Wt{2} = Z2.*V;
       
        % CG: update U
        %-------------
      
        [U,earray1] = CGupdate(b,A, At,D,Dt,Lam,Wt,opts, U, opts.Threshold,opts.CGiterations);
        
        %  Update Lagrange multipliers
        %------------------------------
        
        Lam{1} = Lam{1} - 1.618*opts.beta*(Wt{1} - Dfwd{1});  
        Lam{2} = Lam{2} - 1.618*opts.beta*(Wt{2} - Dfwd{2});    
          
       Dfwd=D(U);
    end 
        %figure(2),imagesc(abs(U)),colormap gray;title(num2str(out));pause(0.1);

        opts.beta=opts.beta*opts.betarate;
       
    end
    
end