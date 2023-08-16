function [freq, pow,time] = BS_MUSIC(Y, beam_w, K, Angle_H,Angle_L)
degrad = pi/180;
[M, N] = size(Y);
Rhat = beam_w*Y * Y'*beam_w' / N;
start = tic;
[U S V] = svd(Rhat);
En = U(:,K+1:end);
GG=En*En';
index = 1;

region = 1;
power = zeros(1,501);
for ag = Angle_L+region:0.01:Angle_H-region
    f = 1/2*sin(ag*degrad);
    A_t = beam_w*exp(1i*2*pi*kron((0:M-1)',f'));
    power(index) =log(1/abs(A_t'*En*En'*A_t));
    index = index +1;
end
power = power/max(power(:));
fi = Peaksearch_m(power, [Angle_L+region:0.01:Angle_H-region]', K);
fi=sort(fi);
% postprocessing and parameter estimation
freq = fi;
pow = 0;
time = toc(start);
end

