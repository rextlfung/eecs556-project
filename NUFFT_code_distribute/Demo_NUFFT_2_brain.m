% Demo to demonstrate the benefit of Optimal least square NUFFT schemes 
%
% [1] M. Jacob,"Optimized least square non uniform fast Fourier transform (OLS-NUFFT)" , 
%     IEEE Transactions of Signal Processing, vol. 57, issue 6, pp. 2165-2177, Feb 2009 
% [2] Z. Yang and M. Jacob,"Mean square optimal NUFFT approximation for efficient non-Cartesian 
%     MRI reconstruction, JMRI, submitted
%  
%  This demo program sets up and compares the different NUFFT schemes in
%  the context of recovering MRI data from 4 channel spiral samples. 
%-------------------------------------------+
% Author: 
% Zhili Yang @ University of Rochester, Mathews Jacob @ University of Iowa
% Email: zhyang@ece.rochester.edu
%        mathews-jacob@uiowa.edu
%-------------------------------------------+


close all;
clear all;

addpath('./NUFFT') 
addpath('./DFT')
addpath('./Optimization')


%% Raw data and trajectory
% Four channel MRI data acquired with spiral trajectory

load('Data/traj_new.mat')
load('Data/rawdatachannel');

% zoom in image
kloc_centered=kloc;

% Shift the image: multiply by phase term in Fourier domain
shift=(exp(1i*2*pi*(shiftX*real(kloc_centered)+shiftY*imag(kloc_centered))));

%% Compute and set up NUFFT forward and backward models: one time task
% These forward models are saved to file to reduce repeated computation
%----------------------------------------------------------------------------

% NUFFT parameters
N=192;
[K Ofactor J ]=deal(N+2,151,4);
diff=K-N;
Order=2;
dcf=ones(1,length(kloc_centered));

path = ['Precomputed_Kernels/',num2str(K),'_',num2str(N),'_',num2str(Ofactor),'_',num2str(Order),'_',num2str(J),'_'];

% Worst case OLS
%----------------

if(~exist([path,'WOLS.mat']))
    [prefilter_WOLS,fn_WOLS]=NUFFT_symmetric(J,N,Ofactor,K,Order);
    save([path,'WOLS.mat'], 'prefilter_WOLS','fn_WOLS');
else
    load([[path,'WOLS.mat']]);
end

A_WOLS = @ (z) NUFFT2_Symmetric(z,K,diff,fn_WOLS,kloc_centered,J,N,Ofactor,prefilter_WOLS,shift);
At_WOLS = @ (z) INUFFT2_Symmetric(z,fn_WOLS,kloc_centered,J,K,N,Ofactor,ones(size(kloc_centered)),prefilter_WOLS,shift);

% MOLS-U
%----------------

if(~exist([path,'MOLS.mat']))
    Energy_distribution = ones(1,N)';
    [prefilter_MOLS,fn_MOLS,error]=NUFFT_Assymetric(J,N,Ofactor,K,Order,Energy_distribution);                                   
    save([path,'MOLS.mat'], 'prefilter_MOLS','fn_MOLS');
else
    load([[path,'MOLS.mat']]);
end

A_MOLS = @ (z)NUFFT2D_general(z,K,fn_MOLS,fn_MOLS,kloc_centered,J,N,Ofactor,prefilter_MOLS,shift);
At_MOLS = @ (z)INUFFT2D_general(z,fn_MOLS,fn_MOLS,kloc_centered,J,K,N,Ofactor,dcf,prefilter_MOLS,shift);


% LS-KB OLS NUFFT
%-----------------------

diff=K-N;
[prefilter_ls,fn_ls] = giveLSInterpolator(N,K,Ofactor,J);
%prefilter_ls_2D=prefilter_ls'*prefilter_ls;
A_LS = @ (z) NUFFT2_Symmetric(z,K,diff,fn_ls,kloc_centered,J,N,Ofactor,prefilter_ls,shift);
At_LS = @ (z) INUFFT2_Symmetric(z,fn_ls,kloc_centered,J,K,N,Ofactor,ones(size(kloc_centered)),prefilter_ls,shift);

% LS-KB with K=2N OLS NUFFT
%-----------------------
Kbig=256;
diffBig=Kbig-N;
Jbig=6;
[prefilter_ls_2N,fn_ls_2N] = giveLSInterpolator(N,Kbig,Ofactor,Jbig);
A_LS_2N = @ (z) NUFFT2_Symmetric(z,Kbig,diffBig,fn_ls_2N,kloc_centered,Jbig,N,Ofactor,prefilter_ls_2N,shift);
At_LS_2N = @ (z) INUFFT2_Symmetric(z,fn_ls_2N,kloc_centered,Jbig,Kbig,N,Ofactor,ones(size(kloc_centered)),prefilter_ls_2N,shift);



%% Solve  inverse problems


% Optimization parameters
%------------------------

opts.mu = 1e-7;            % regularization parameter
opts.beta=1;               % Initial beta value; scalar used for penalty term
opts.betarate = 1.5;         % Increment use for beta
opts.ContThreshold = 1e-2; % Threshold for updating beta; increment beta when error saturates
opts.iter = 5;            % no of iterations
opts.outiter = 10;
opts.CGiterations = 45;    % Number of maximum CG iterations 
opts.Threshold = 1e-9;     % Exit condition

close all;
pause(1);

%
% Worst case Optimum Least Square NUFFT
%---------------------------------------
   
    % Reconstruct all 4 channels independently
    for ch=1:4
        y=At_WOLS(transpose(X{ch}));
        [x1(:,:,ch),cost,datacons, tvnorm] = TVrecon(A_WOLS,At_WOLS,y,transpose(X{ch}),opts);
    end
    
    % Sum of squares reconstruction
    SOS_WOLS=sqrt(sum(abs(x1).^2,3));
    figure(1);imshow(flipud(abs(SOS_WOLS')),[],'InitialMagnification','fit');axis off
    title(['WOLS']);
    pause(1);


% Mean sqare Optimum Least Square NUFFT assuming uniform energy distribution
%---------------------------------------------------------------------------

    % Reconstruct all 4 channels independently
    for ch=1:4
        y=At_MOLS(transpose(X{ch}));
        [x_u(:,:,ch),cost,datacons, tvnorm] = TVrecon(A_MOLS,At_MOLS,y,transpose(X{ch}),opts);
    end
    
    % Sum of squares reconstruction
    SOS_MOLS=sqrt(sum(abs(x_u).^2,3));
    figure(2);imshow(flipud(abs(SOS_MOLS')),[],'InitialMagnification','fit');axis off
    title(['MOLS']);
    pause(1);

% LS-KB OLS NUFFT
%-----------------

    
    % Reconstruct all 4 channels independently
    for ch=1:4
        y=At_LS(transpose(X{ch}));
        [x_ls(:,:,ch),cost,datacons, tvnorm] = TVrecon(A_LS,At_LS,y,transpose(X{ch}),opts);
    end
    % Sum of squares reconstruction

    SOS_KB=sqrt(sum(abs(x_ls).^2,3));
    figure(3);imshow(flipud(abs(SOS_KB')),[],'InitialMagnification','fit');axis off
    title(['LS-KB']);
pause(1);


% LS-KB OLS NUFFT K=2N
%-----------------

    
    % Reconstruct all 4 channels independently
    for ch=1:4
        y=At_LS_2N(transpose(X{ch}));
        [x_ls_2N(:,:,ch),cost,datacons, tvnorm] = TVrecon(A_LS_2N,At_LS_2N,y,transpose(X{ch}),opts);
    end
    % Sum of squares reconstruction

    SOS_KB_2N=sqrt(sum(abs(x_ls_2N).^2,3));
    figure(4);imshow(flipud(abs(SOS_KB_2N')),[],'InitialMagnification','fit');axis off
    title(['LS-KB']);
pause(1);