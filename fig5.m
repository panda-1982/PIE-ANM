% test of SPA in the ULA case

clear all
% close all
% clc


K = 3;          % source number
M = 65;         % array length
B = 6;          % beamspace dimension
N = 50;        % snapshot number
sigma = 1;      % noise variance
rmes = 0;
degrad = pi/180;
% true frequency and power
theta = [15;18.5;22];
theta_i = [46];
K_I = 1;
f = 1/2*sin(theta*degrad);
f_i = 1/2*sin(theta_i*degrad);
p = ones(K+K_I,1);
rmse_ie = 0;
rmse_ie1 = 0;
rmse_ie_bm = 0;
rmse_ie_bm1 = 0;
rmse_FS = 0;
rmse_wane = 0;
rmse_music = 0;
rmse_gl = 0;
snr_index = 1;
mse_ie = zeros(1,11);
mse_FS = zeros(1,11);
mse_ie1 = zeros(1,11);
mse_wane = zeros(1,11);
mse_ie_bm = zeros(1,11);
mse_ie_bm1 = zeros(1,11);
mse_music = zeros(1,11);
mse_gl = zeros(1,11);
time_cost_ie = 0;
time_cost_ie1 = 0;
time_cost_FS = 0;
time_cost_gl = 0;
time_cost_music = 0;
time_cost_ie_bm = 0;
time_cost_wane = 0;
time_cost_ie_bm1 = 0;
alpha = 0.5;
A_e = exp(1i*2*pi*kron((0:M-1)',f'));
A_i = exp(1i*2*pi*kron((0:M-1)',f_i'));
S = sqrt(diag(p))*exp(1i*2*pi*rand(K+K_I,N));
S(2,:) = alpha*S(1,:)+(1-alpha)*S(2,:);
S(3,:) = exp(1i*rand(1)*2*pi)*alpha*S(1,:)+(1-alpha)*S(3,:);


iter = 300;
Beam_w = (1/M)*exp(1i*2*pi/M*(-13:B-14)'*(-(M-1)/2:(M-1)/2));


ag_grid = asin(linspace(-1,1,3601))/degrad;
A_p = zeros(1,3601);
for index = 1:3601
    f = 1/2*sin(ag_grid(index)*degrad);
    A_t = Beam_w*exp(1i*2*pi*kron((-(M-1)/2:(M-1)/2)',f'));
    A_p(index) = sqrt(abs(A_t'*A_t));
end
max_p = sqrt(max(A_p));
A_p = A_p/(max_p^2);

Angle_H = 25.5;
Angle_L = 12.5;

angle_index = A_p>0.7;
angle_regin = ag_grid(angle_index);
Angle_H1 = max(angle_regin);
Angle_L1 = min(angle_regin);

ag_grid = linspace(sin(Angle_L*degrad),sin(Angle_H*degrad),3601);
ag_mainlobe = (ag_grid>sin(Angle_L1*degrad)).*( ag_grid<sin(Angle_H1*degrad));
index_H1 = find(ag_mainlobe,1,'last');
index_L1 = find(ag_mainlobe,1,'first');


for index = 1:3601
    f = 1/2*ag_grid(index);
    A_t = Beam_w*exp(1i*2*pi*kron((0:M-1)',f'));
    A(:,index) = A_t;
    index = index + 1;
end

g = zeros(B,3601);
g = zeros(B,3601);
for num = 1:B
    g(num,:)=fftshift(ifft(A(num,:),3601));
end
index = -1800:1800;
g_sum = sum(abs(g),1);
g_max = max(g_sum);
index_pick = index(g_sum>0.05*g_max);
bias = ceil((max(index_pick)-min(index_pick))/2)+min(index_pick);

M_order = 10;
ag1_grid = linspace(-pi,pi,3601)-pi;
ag2_grid = (0:3600)/3601*2*pi;
G=g(:,1801-M_order:1801+M_order);

ag_mainlobe = ag_grid(index_L1:index_H1);
f_mainlobe = 1/2*ag_mainlobe;
d_f = 1/2*degrad*cos(asin(ag_mainlobe));
A_t_mainlobe = Beam_w*exp(1i*2*pi*kron((-(M-1)/2:(M-1)/2)',f_mainlobe));

ag_sidelobe = [ag_grid(1:index_L1-1) ag_grid(index_H1+1:end)];
f_sidelobe = 1/2*ag_sidelobe;
A_t_sidelobe = Beam_w*exp(1i*2*pi*kron((-(M-1)/2:(M-1)/2)',f_sidelobe));

A_t_mainlobe_coef_samp = A_t_mainlobe(:,1:10:end);
A_t_sidelobe_samp = A_t_sidelobe(:,1:10:end);

M_order = 6;
M11 = M_order-bias;
M12 = M_order+bias;
d = exp(-1i*(-M11:M12)'*ag1_grid);
d_sidelobe = d(:,[1:index_L1-1 index_H1+1:end]).*exp(-1i*pi*(M-1)/2*ag_sidelobe);
d_mainlobe = d(:,index_L1:index_H1).*exp(-1i*pi*(M-1)/2*ag_mainlobe);
d_sidelobe_samp = d_sidelobe(:,1:10:end);

% cvx_solver mosek;

cvx_begin
variable G1(B,2*M_order+1) complex,
norm(G1*d_mainlobe-A_t_mainlobe,'fro')<= .2

minimize norm(d_sidelobe_samp'*G1'*A_t_mainlobe_coef_samp,'fro')
cvx_end

A_apr1 = G1*d;




d = 20;
dd = (d-1)/2;
t = 0:M-1;
t_t = t/(M-1);
f_grid = linspace(-1,1,3601);
gamma = ((min(ag_grid))+(max(ag_grid)))/2;
delta = ((max(ag_grid))-(min(ag_grid)))/2;
c = pi*(M-1)*delta;

Psi = zeros(2*dd+1,2*dd+1);
h = zeros(2*dd+1,M);
for num=1:(2*dd+1)
    P = pswf(num-1,c);
    Psi(:,num) = P((-dd:dd)/dd);
%     h(num,:) = P(-t_t);
    h(num,:) = P(2*t_t-1);
end
DD = diag(exp(1i*pi*gamma*(0:M-1)));
v = exp(1i*c*(0:(2*dd))'/(2*dd)*(f_grid));
G_wane = Beam_w*DD*h.';
w_l = c/(2*dd);
J1 = [eye(2*dd) zeros(2*dd,1)];
J2 = [zeros(2*dd,1) eye(2*dd)];
J_add = tan(w_l/2)*(J1+J2);
J_sub = tan(w_l/2)*(J1-J2);
vv = exp(1i*c*(0:(2*dd))'/(2*dd)*(sin(theta'*degrad)-gamma)/delta);
norm(G_wane*inv(Psi)*vv-Beam_w*A_e,'fro')
I = 1/2*sin([-Angle_H,-Angle_L]*degrad);

theta

SNR = linspace(-5,15,11);
load s_fig2;
inr = 15;
max_iter = 2000;
tau = 0.02;
for iS = 1:length(SNR)
    snr = SNR(iS);
%     snr = 0;
    parfor iter_num=1:iter
        warning off;
        p_s = ones(1,K);
        p_i = 10^(inr/20);
        Yt        = [A_e A_i]*diag([p_s p_i])*S;
        Y = awgn(Yt,snr,'measured');
        Rhat = Beam_w*Y*Y'*Beam_w'/N;
        Rhat = B*Rhat/sum(real(diag(Rhat)));
        delta1 = svd(Rhat);
        delta_min = delta1(end);
        delta_max = delta1(1);
        eta =.1*(delta1);
        [freq_music,pow_music,time_music] = BS_MUSIC(Y,Beam_w,K,Angle_H,Angle_L);
        [freq_wane,pow,time_wane] = WANE(Rhat, G_wane, Psi,K, J_add, J_sub,c,delta,gamma);
        [freq_GL,pow_music,time_GL] = GL_ANM(Y,Beam_w,K);   
        freq_GL = sort(-freq_GL);
        [freq_IE, pow,time_IE] = IE_GL_ANM(Rhat, Beam_w, G, Angle_H,Angle_L,K);
        [freq_IE1, pow1,time_IE1] = IE_GL_ANM(Rhat, Beam_w, G1, Angle_H,Angle_L,K);
        [freq_admm_BM,Tu,time_admm_IE] = IE_GL_admm_BM1(Rhat,Beam_w,G,-Angle_L,-Angle_H,K,tau,max_iter);
        [freq_admm_BM1,Tu1,time_admm_IE1] = IE_GL_admm_BM(Rhat,Beam_w,G1,-Angle_L,-Angle_H,K,tau,max_iter);
        freq_admm_BM = sort(-freq_admm_BM);
        freq_admm_BM1 = sort(-freq_admm_BM1);
        [freq_FS, pow,time_FS] = FS_ANM(Rhat, Beam_w, I, K,eta);
        freq_FS = sort(-freq_FS);
        time_cost_ie = time_cost_ie + time_IE(end);
        time_cost_ie1 = time_cost_ie1 + time_IE1(end);
        time_cost_ie_bm = time_cost_ie_bm + time_admm_IE;
        time_cost_ie_bm1 = time_cost_ie_bm1 + time_admm_IE1;
        time_cost_FS = time_cost_FS + time_FS(end);
        time_cost_gl = time_cost_gl + time_GL(end);
        time_cost_wane = time_cost_wane + time_wane(end);
        time_cost_music = time_cost_music + time_music;
        
        % MSE
        rmse_music = rmse_music+sum((theta - sort(freq_music(1:K))).^2);
        rmse_ie = rmse_ie+sum((theta - sort(freq_IE(1:K))).^2) ;
        rmse_ie1 = rmse_ie1+sum((theta - sort(freq_IE1(1:K))).^2) ;
        rmse_FS = rmse_FS+sum((theta - sort(freq_FS(1:K))).^2) ;
        rmse_gl = rmse_gl +sum((theta - sort(freq_GL(1:K))).^2) ;
        rmse_ie_bm = rmse_ie_bm+sum((theta - sort(freq_admm_BM(1:K))).^2);
        rmse_ie_bm1 = rmse_ie_bm1+sum((theta - sort(freq_admm_BM1(1:K))).^2); 
        rmse_wane = rmse_wane+sum((theta - sort(freq_wane(1:K))).^2); 

    end
    mse_ie(snr_index)=sqrt(rmse_ie/(iter*K));
    mse_ie1(snr_index)=sqrt(rmse_ie1/(iter*K));
    mse_FS(snr_index)=sqrt(rmse_FS/(iter*K));
    mse_music(snr_index)=sqrt(rmse_music/(iter*K));
    mse_gl(snr_index)=sqrt(rmse_gl/(iter*K));
    mse_wane(snr_index)=sqrt(rmse_wane/(iter*K));
    mse_ie_bm(snr_index) = sqrt(rmse_ie_bm/(iter*K));
    mse_ie_bm1(snr_index) = sqrt(rmse_ie_bm1/(iter*K));
%     sqrt(rmse_wane/(iter*K))
%     sqrt(rmse_gl/(iter*K))
%     sqrt(rmse_ie/(iter*K))
%     sqrt(rmse_ie1/(iter*K))
%     sqrt(rmse_ie_bm/(iter*K))
%     sqrt(rmse_ie_bm1/(iter*K))
%     sqrt(rmse_music/(iter*K))
%     time_cost_ie1
%     time_cost_ie
%     time_cost_ie_bm
%     time_cost_ie_bm1
%     time_cost_wane
%     time_cost_music
    
    snr_index = snr_index + 1;
    rmse_ie = 0;
    rmse_ie1 = 0;
    rmse_music = 0;
    rmse_gl = 0;
    rmse_ie_bm = 0;
    rmse_ie_bm1 = 0;
    rmse_wane = 0;
    rmse_FS = 0;
end
snr = SNR;
mz = 6;
lw = 0.5;

h=figure;
semilogy(snr,mse_music, '-p', 'markersize', mz, 'linewidth', 2); hold on;
 
semilogy(snr,mse_gl, '-d', 'markersize', mz, 'linewidth', 2);
semilogy(snr,mse_FS, '-*', 'markersize', mz, 'linewidth', 2);
semilogy(snr,mse_wane, '-x', 'markersize', mz, 'linewidth', 2);
semilogy(snr,mse_ie, '-o', 'markersize', mz, 'linewidth',2);
semilogy(snr,mse_ie1, '--^', 'markersize', mz, 'linewidth', 2);
semilogy(snr,mse_ie_bm, '->', 'markersize', mz, 'linewidth', 2);
semilogy(snr,mse_ie_bm1, '-.s', 'markersize', mz, 'linewidth', 2);
grid on;
legend('MUSIC','ANM','FS-ANM','WANE','PIE-ANM (25)','PIE-ANM (27)','PIE-BM-ADMM (25)','PIE-BM-ADMM (27)')
xlabel('SNR (dB)')
ylabel('RMSE(\circ)')
axis([-5 15 0.04 2]);
savefig(h,'SNR.fig')
close(h)

h=figure;
time = [time_cost_music time_cost_ie_bm1 time_cost_ie_bm  time_cost_ie1 time_cost_ie time_cost_wane time_cost_gl time_cost_FS]/(11*iter);
c = categorical({'Figure. 5','Figure. 5_{repeat}'});
bar(c,[time;time]);
legend('MUSIC','PIE-BM-ADMM (25)','PIE-BM-ADMM (27)','PIE-ANM (25)','PIE-ANM (27)','WANE','ANM','FS-ANM')
ylabel('CPU Time(s)')
grid on;
savefig(h,'Time.fig')
close(h)