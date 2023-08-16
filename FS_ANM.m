function [freq,pow,time] = FS_ANM(yonOmega, Omega, I, K,eta)

% [q, y, f_root, f_VD, p_VD] = FS_ANM(yonOmega, Omega, I, N, eta)
%
% FS_ANM implements the frequency-selective atomic norm minimization (FS-ANM).
%
% Inputs:
%   yonOmega: samples indexed by Omega
%   Omega: index set
%   I: union of two intervals on the unit circle, J*2 matrix
%   N: length of full data
%   eta: upper bound on the Frobenius norm of noise
% Output:
%   q: coefficients of the dual polynomial
%   y: recovered full data
%   f_root: recovered frequencies computed using the root-finding method
%   f_VD: recovered frequencies computed using the Vandermonde
%   decomposition method
%   p_VD: recovered magnitude
%
% Written by Zai Yang, April 2016
% References:
% [1] Z. Yang and L. Xie, "Vandermonde decomposition on intervals and its use for continuous compressed sensing," 35th Chinese Control Conference (CCC), pp. 4877-4882, July 2016.
% [2] Z. Yang and L. Xie, "Frequency-selective Vandermonde decomposition of Toeplitz matrices with applications," arXiv: 1605.02431, 2016.

if nargin < 5
    eta = 0;
end
if nargin < 4
    N = Omega(end);
end

% compute r0 & r1 based on I
r0 = -2*cos(pi*(I(2) - I(1))) * sign(I(2) - I(1));
r1 = exp(-1i*pi*sum(I)) * sign(I(2) - I(1));

[B , N] = size(Omega);

% eta = 0.01;
cvx_tic;
cvx_quiet true
cvx_precision default
% cvx_solver mosek
cvx_begin sdp

variable y(N,B) complex,
variable x(B,B) hermitian,
variable t(N,1) complex,
dual variable X


norm(Omega*y - yonOmega,'fro') <= eta;

[x, y'; y, toeplitz(t)] >= 0 : X,
r1*toeplitz([t(2); conj(t(1:N-2))], t(2:N)) + r0*toeplitz(t(1:N-1))...
    + conj(r1)*toeplitz(conj(t(2:N)), [conj(t(2)); t(1:N-2)]) >= 0,

minimize trace(x) + real(trace(toeplitz(t)));

cvx_end

% frequency retrieval using Vandermonde decomposition
T = toeplitz(t);
% [f_VD, p_VD] = FS_VD(T, I);
degrad = pi/180;
[U SS V] = svd(T);
En = U(:,K+1:end);
GG=En*En';
MM=N+1;
a = zeros(2*MM-1,1);
for j=-(MM-1):(MM+1)
    a(j+MM) = sum( diag(GG,j) );
end





ra=roots([a]);
rb=ra(abs(ra)<1);

rc = rb;
% pick the n roots that are closest to the unit circle
[~,I]=sort(abs(abs(rc)-1));
if size(I,1)>=K
    w=angle(rc(I(1:K)));
else
    w=angle([rc ;zeros(K-size(I,2),1)]);
end
freq=sort(asin(w/pi)/degrad);
time = cvx_toc;
pow =  T;
end

