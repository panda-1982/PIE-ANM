function [Tu,error] = admm_BM(Y,B,G,K,tau,rho_1,rho_2,r,maxIter)
%size of mesurements
Z = Y;
Z_size = size(Z);
m = Z_size(1,1);
n = size(G,2);
L = Z_size(1,2);

%parameters
tol_abs = 2e-4; %absolute tolerance 1e-3
tol_rel = 2e-4; % iteration level relative tolerance 1e-4
converged = 0;

%initialization
YOld = zeros(n+L)+1i*zeros(n+L);
Lambda = zeros(n+L)+1i*zeros(n+L);
Q = zeros(n+L,r)+1i*zeros(n+L,r);
Q_hat=Q;
Lambda1 = rand(n+L,r);

normalizer = 1./[n; ((n-1):-1:1)'];
e1 = zeros(n,1); e1(1)=1;
Beam=B;
B = G;
B_rho = inv(B'*B+2*rho_1*eye(n));
error = zeros(1,maxIter);
for count = 1:1:maxIter
    
    %update the variables W,X, and u
    W = YOld(n+1:n+L,n+1:n+L)...
        + (1/rho_1)*(Lambda(n+1:n+L,n+1:n+L) - (tau/2)*eye(L));
    X = (B_rho)*(B'*Z+ 2*Lambda(n+1:n+L,1:n)' + rho_1*YOld(1:n,n+1:n+L)+ rho_1*YOld(n+1:n+L,1:n)');
    u = (1/rho_1) * normalizer.*conj(toeplitz_adjoint(Lambda(1:n,1:n))...
        + rho_1*toeplitz_adjoint(YOld(1:n,1:n)) - (tau/2)*n*e1);
    
    temp = [toeplitz(conj(u)), X; X', W];
    
    
    
    
    Q = (rho_2*Q_hat+rho_1*temp*Q_hat-Lambda*Q_hat-Lambda1)*inv(rho_1*Q_hat'*Q_hat+rho_2*eye(r));
    Q_hat = (rho_2*Q+rho_1*temp'*Q-Lambda'*Q+Lambda1)*inv(rho_1*Q'*Q+rho_2*eye(r));
    QQ = Q*Q_hat';
    Y = (QQ+QQ')/2;
    
    %stop criteria.
    pri_res = temp - Y;
    
    dual_res_norm = rho_1*sqrt((norm(Y(1:n,n+1:n+L)-YOld(1:n,n+1:n+L),'fro'))^2/L+(norm(toeplitz_adjoint(Y(1:n,1:n)-YOld(1:n,1:n)),'fro'))^2+(norm(Y(n+1:n+L,n+1:n+L)-YOld(n+1:n+L,n+1:n+L),'fro'))^2/L);
    dual_var_adj_norm = rho_1*sqrt((norm((Lambda(1:n,n+1:n+L)+Lambda(n+1:n+L,1:n)')/2,'fro'))^2/L+(norm(toeplitz_adjoint(Lambda(1:n,1:n)),'fro'))^2+(norm(Lambda(n+1:n+L,n+1:n+L),'fro'))^2/L);
    
    pri_tol = (n+L)*tol_abs + tol_rel*max(norm(temp,'fro'),norm(Y,'fro'));
    dual_tol = sqrt(2*n+1)*tol_abs + tol_rel*dual_var_adj_norm;
    
    err_rec_primal = norm(pri_res,'fro');
    err_rec_dual = dual_res_norm;
    
    converged = and(err_rec_primal<pri_tol,err_rec_dual<dual_tol);
    Tu = toeplitz(conj(u));
    error(count) = tau/2*real(trace(Tu)+trace(W))+1/2*norm(B*X-Z,'fro');
    if converged, break; end
    Lambda = Lambda + rho_1*(Y - temp);
    Lambda1 = Lambda1+rho_2*(Q-Q_hat);
    YOld = Y;
end
Tu = toeplitz(conj(u));



function T = toeplitz_adjoint(A)
N = size(A,1);
T = zeros(N,1)+1i*zeros(N,1);
T(1) = sum(diag(A));
for n = 1:(N-1)
    T(n+1) = sum(diag(A,n));
end


   
