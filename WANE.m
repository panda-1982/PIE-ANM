function [freq,pow,time] = WANE(Rhat, G_wane, Psi,K, J_add, J_sub,c,delta,gamma)






[B L]=size(G_wane);
    Rhat = B*Rhat/sum(real(diag(Rhat)));
    delta1 = svd(Rhat);
    delta1 = delta1(end);
    eta =.1*(delta1);
dd = (L-1)/2;
degrad = pi/180;
cvx_tic;
cvx_quiet true;
cvx_precision default
% cvx_solver mosek
cvx_begin sdp
  
  variable y(L,B) complex,
  variable x(B,B) hermitian,
  variable t(L,1) complex,
  dual variable X
  

  norm(G_wane*y - Rhat,'fro') <= eta;

  [x, y'*Psi'; Psi*y, toeplitz(t)] >= 0 : X,
  J_add* toeplitz(t)*J_add'-J_sub* toeplitz(t)*J_sub' >= 0,
  
 minimize trace(x) + real(trace(toeplitz(t)));
  
cvx_end
[U SS V] = svd(toeplitz(t));
En = U(:,K+1:end);
GG=En*En';
MM=L+1;
a = zeros(2*MM-1,1);
for j=-(MM-1):(MM+1)
    a(j+MM) = sum( diag(GG,j) );
end

ra=roots([a]);
rb=ra(abs(ra)<1);
[~,I]=sort(abs(abs(rb)-1));
w=angle(rb(I(1:K)));
freq = sort(-asin(w/c*(2*dd)*delta-gamma)/degrad);
    time = cvx_toc;
    pow =  toeplitz(t); 
end    

