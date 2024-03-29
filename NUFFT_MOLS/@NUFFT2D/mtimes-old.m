function Ax = mtimes(A, x)
% spiral_nufft_matrix_MC_GPU/MTIMES   
% Implement Ax for Polynomial&Euclidian decomposition 
% first then b0 field and multiplication by sensitivity map for each
% channel then fft for 3-D data x of size N1xN2x(N3+Np)
% before nufft, Ax has size N1xN2xN3xNc with Nc as the no. of channels

% x = [xpeak xpol] where xpeak is the baseline-removed x for only frange; 
% size(mask)=size(x) xpol is the coeffs of x decomposed in polynomial 
% space along f (or n3) dimension
% adapted from spiral_nufft_matrix_GPU/mtimes

[N1, N2, N3] = size(A.mask); 
%for reshape the parameters should be in CPU
N1 = single(N1); N2 = single(N2); N3 = double(N3);
Nc = size(A.sense, 3);

    if A.vectorize
        x = reshape(x, N1, N2, N3);
    end
    x = x.*single(A.mask);
    x = (reshape(x, N1*N2, N3)).';
    
    % Transforming to time domain
    %x = ifft(x, [], 3); 
    %x = ifft(x);

    % Shifting the data
    x = reshape(x.', N1, N2, N3);
    x = x.*A.w;
    
    %x = b0T2comp3D(x, A.w, A.T);

    % Accounting for coil sensitivities and spiral trajectory
    
    AxTemp = zeros(N1,N2,N3,Nc,'single');
    for(ch = 1:Nc)
        temp = zeros(2*N1,2*N2,N3);
        temp(N1/2+1:end-N1/2,N2/2+1:end-N2/2,:) = Sens3D(x, A.sense(:,:,ch));
        temp = ifftshift(ifft2(fft2(fftshift(temp)).*A.Weight));
        temp = Sens3D(temp(N1/2+1:end-N1/2,N2/2+1:end-N2/2,:),conj(A.sense(:,:,ch)));
        AxTemp(:,:,:,ch) = temp;
        %NUFFT3_matrix(Ax1, A.S, A.Nuparam(1), A.Nuparam(2), A.prefilter_2D); %3D NUFFT 
    end
    Ax = zeros(N1,N2,N3, 'single');
    for ch = 1:Nc
        Ax = Ax + squeeze(AxTemp(:,:,:,ch));
    end
    clear AxTemp;
    Ax = Ax.*conj(A.w);
    %Ax = b0T2comp3D(Ax, -A.w, A.T);
    Ax = (reshape(Ax, N1*N2, N3)).';
    %Ax = fft(Ax)/single(N3);
    Ax = reshape(Ax.', N1, N2, N3);%fft(x,[],3)/(N3-Np);
    
   
    Ax = Ax.*single(A.mask); % this is required
    if A.vectorize
        Ax = Ax(:); % suppose that M1=N1, M2=N2
    end


function raw = INUFFT3_matrix(Raw, S, K, N, prefilter_2D)
K_G = single(K);
N3 = size(Raw,2);
S = ctranspose(S);

pdata = zeros(K,K,N3, 'double');
%pdata = repmat(reshape(Raw(:,1),K,K),[1,1,N3]);

Raw = double(Raw);
for n = 1:N3
    data = Raw(:,n);
    data = S*data(:);%real(data(:)) + (S*imag(data(:)))*i;
    pdata(:,:,n) = reshape(data,K,K);
end

pdata = single(pdata);

pdata = ifftshift(ifft2(ifftshift(pdata)))*K_G;

pdata = pdata.*repmat(prefilter_2D,[1 1 N3]);
% raw is in (kx,ky,n3) and Raw in (data, n3)

diff = K-N;
raw = pdata(diff/2+1:end-diff/2, diff/2+1:end-diff/2, :); 


function Raw = NUFFT3_matrix(raw, S, K, N, prefilter_2D)

N3 = size(raw,3);

paddedimage = zeros(K, K, N3, 'single');
diff = K-N; %diff should be at cpu
paddedimage(diff/2+1:end-diff/2,diff/2+1:end-diff/2,:) = raw;

paddedimage = paddedimage.*repmat(prefilter_2D,[1 1 N3]);

KFFTimag = fftshift(fft2(fftshift(paddedimage)));
clear paddedimage
Raw = zeros(size(S,1), size(raw,3), 'double');
 
KFFTimag = double(KFFTimag);
for n = 1:N3 % cpu version is faster since S is sparse
    KF = KFFTimag(:,:,n);
    Raw(:,n) = S*KF(:)/K; %(S*real(KF(:)) + (S*imag(KF(:)))*i)/K_G;
end

% B0 inhomogeniety compensation 
function x = b0T2comp3D(x, w0, T2)

%x = ifft(x, [], 3);
%for k = 2:size(x, 3)
%    x(:,:,k) = x(:,:,k).*exp((T2(:,:) + j*w0(:,:))*(k-1));
%end
[N1 N2 N3] = size(x);
x = x(:);
Tw = T2(:) + j*w0(:);
N = N1*N2;
% Tw = repmat(Tw,[N3 1]); %this gives erroneous result on gpu since
% size(Tw)>1MB this is a repmat bug in jacket
% we do this other way similar to repmat
% function (we also can do this for K below but not required since each 
% size of K (rows or columns) is less than 1MB)
mind = single(1:N)';
mind = mind(:, ones(1,N3));
Tw = Tw(mind,1);

K = repmat(single(0:N3-1), [N,1]);
Tw = Tw.*K(:);

x = x.*exp(Tw); %this is much faster for GPU maybe 50times
x = reshape(x, N1, N2, N3); 


% sensitivity compensation 
function x = Sens3D(x, sens)

x = x.*repmat(sens, [1, 1, size(x,3)]);

