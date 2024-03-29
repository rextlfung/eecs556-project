%--------------------------------------------------------------------------
% CG solution to region limited problem
% X,earray1] = CGupdate(b,A, At,D,Dt,Lam,Y,C, X, THRESHOLD,Niter)
%--------------------------------------------------------------------------

function [X,earray1] = CGupdate(b,A, At,D,Dt,Lam,Y,C, X, THRESHOLD,Niter)

oldcost = 0;
earray1 = [];
lam = 0.5*C.mu*C.beta;

LamTV = Lam; 
eTV=double(0);
for i=1:Niter,
    
    resY = (A(X) - b);
    eY = sum(abs(resY(:)).^2);
    
    Dx = D(X);
    resTV{1} = Dx{1}-Y{1};
    resTV{2} = Dx{2}-Y{2};
    
    LamTV{1}=LamTV{1}+1i*1e-18;
    LamTV{2}=LamTV{2}+1i*1e-18;
    
    eTV = lam*(abs(sum(conj(resTV{1}(:)).*resTV{1}(:))) + abs(sum(conj(resTV{2}(:)).*resTV{2}(:))) );
    eTV = eTV + C.mu *abs(conj(sum(conj(LamTV{1}(:).*resTV{1}(:))))+conj(sum(conj(LamTV{2}(:).*resTV{2}(:)))));
    
    
    cost1 = eY + eTV;
    
    earray1 = [earray1,cost1];
    
    if(abs(cost1-oldcost)/abs(cost1) < THRESHOLD)
        %display('exiting')
        break;
    end
    oldcost = cost1;
    
  %  conjugate gradient direction
   % ------------------------------
    
    % gradient: gn
    
    gn = At(A(X)-b) + lam*Dt(resTV);
    gn = 2*gn+C.mu*Dt(LamTV);
    
    % search direction: sn  
    if(i==1)
        sn = gn;                                          
        oldgn = gn;
    else
        gamma = abs(sum(sum(conj(gn(:)).*gn(:)))/sum(sum(conj(oldgn(:)).*oldgn(:))));
        sn = gn + gamma*sn; 
        oldgn = gn;
    end
    
    % line search
    %-------------
    Asn = A(sn);  
    Dsn = D(sn);
    
    numer = sum(conj(Asn(:)).*resY(:)) + lam*sum(sum(sum(conj(resTV{1}).*Dsn{1})))  + lam*sum(sum(sum(conj(resTV{2}).*Dsn{2})));
    numer = numer +  0.5*C.mu*sum(sum(sum(conj(LamTV{1}).*Dsn{1}))) + 0.5*C.mu*sum(sum(sum(conj(LamTV{2}).*Dsn{2}))); 
    denom = sum(conj(Asn(:)).*Asn(:)) + lam*sum(sum(sum(conj(Dsn{1}).*Dsn{1}))) + lam*sum(sum(sum(conj(Dsn{2}).*Dsn{2}))); 

    if(denom < 1e-18)
                %display('denom exiting')

        break;
    end
    alpha = -real(numer)/real(denom);
   
    % updating
    %-------------
    
    X = (X + alpha*sn);
end

    
