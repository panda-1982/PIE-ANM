function [freq, pow,time] = GL_ANM(Y, beam_w, K)

% Written by Pan Jie, 2019   E-mail:panjie@yzu.edu.cn


[M, N] = size(Y);



    Rhat = Y * Y' / N;

cvx_quiet true
cvx_precision default

    [freq, pow,time] = SPA_nonsing(Rhat, beam_w, N,K);
   
end



function [freq, pow,time] = SPA_nonsing(Rhat,beam_w, N, K)
% SPA when Rhat is nonsingular
    M = size(Rhat, 1);
    B = size(beam_w,1);
    Rhat = beam_w*Rhat*beam_w';
    Rhat = B*Rhat/sum(real(diag(Rhat)));
    delta = svd(Rhat);
    delta = delta(end);
    eta =.3*(delta);
    degrad = pi/180;
    
    
    
%     cvx_solver sdpt3
%  cvx_solver mosek
%     cvx_quiet false;

    cvx_tic
    cvx_begin sdp
    variable x(B,B) complex,
    variable u(M) complex,
    variable R(M,B) complex
    
    [x R'; R toeplitz(u)] >= 0,
    norm(beam_w*R-Rhat,'fro')<=eta
    minimize real(trace(x)) + real(trace(toeplitz(u)));
    cvx_end
    sig = 0;
    R_cov = toeplitz(u);
    [U SS V] = svd(R_cov);
%     svd(beam_w*R_cov*beam_w')
    En = U(:,K+1:end);
    GG=En*En';
    MM=M+1;
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
    pow = 0;
    time = cvx_toc;
end







