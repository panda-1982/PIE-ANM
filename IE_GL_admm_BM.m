function [freq,Tu,time] = IE_GL_admm_BM(Y,B,G,Angle_H,Angle_L,K,tau,max_iter)

%size of mesurements
start=tic;

r=3*K;
[Tu] = admm_BM(Y,B,G,K,tau,.1,2,r,max_iter);
degrad = pi/180;
region = 0;
R_cov = G*Tu*G';
[U, ~, ~] = svd(R_cov);
En = U(:,K+1:end);
GG=B'*(En*En')*B;
NN = size(B,2);
MM=NN+1;
a = zeros(2*MM-1,1)+1i*zeros(2*MM-1,1);
for j=-(MM-1):(MM+1)
    a(j+MM) = sum( diag(GG,j) );
end

ra=roots([a]);
rb=ra(abs(ra)<1);
freq_bsanm=asin(angle(rb)/pi)/degrad;
rc_i = freq_bsanm<(Angle_H-region)&freq_bsanm>(Angle_L+region);
rc = rb(rc_i);
% pick the n roots that are closest to the unit circle
[~,I]=sort(abs(abs(rc)-1));
if size(I,1)>=K
    w=angle(rc(I(1:K)));
else
    w=angle([rc ;zeros(K-size(I,2),1)]);
end
freq=sort(asin(w/pi)/degrad);
time = toc(start);

