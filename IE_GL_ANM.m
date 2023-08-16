function [freq, pow,time] = IE_GL_ANM(Rhat, beam_w, W, Angle_H,Angle_L,K)
%
% Written by Pan Jie, 2022   E-mail:panjie@yzu.edu.cn


[B, M] = size(W);
NN = size(beam_w,2);
M = (M-1)/2;
Rhat = B*Rhat/sum(real(diag(Rhat)));
delta = svd(Rhat);
delta = delta(end);
eta =.1*(delta);
% eta = 0.01;
degrad = pi/180;
% cvx_solver mosek
cvx_tic
cvx_quiet true
cvx_precision default

cvx_begin sdp quiet
variable x(B,B) hermitian,
variable u(2*M+1) complex,
variable R(2*M+1,B) complex

[x R'; R toeplitz(u)] >= 0,
norm(W*R-Rhat,'fro')<=eta
minimize trace(x) + real(trace(toeplitz(u)));
cvx_end

region = 1;


R_cov = W*toeplitz(u)*W';
[U SS V] = svd(R_cov);
En = U(:,K+1:end);
GG=beam_w'*En*En'*beam_w;
MM=NN+1;
a = zeros(2*MM-1,1);
for j=-(MM-1):(MM+1)
    a(j+MM) = sum( diag(GG,j) );
end

ra=roots([a]);
rb=ra(abs(ra)<1);
freq_bsanm=asin(angle(rb)/pi)/degrad;
rc_i = freq_bsanm>(-Angle_H+region)&freq_bsanm<(-Angle_L-region);
rc = rb(rc_i);
% pick the n roots that are closest to the unit circle
[~,I]=sort(abs(abs(rc)-1));
if size(I,1)>=K
    w=angle(rc(I(1:K)));
else
    w=angle([rc ;zeros(K-size(I,2),1)]);
end
freq=sort(-asin(w/pi)/degrad);
time = cvx_toc;
pow =  R_cov;
end











