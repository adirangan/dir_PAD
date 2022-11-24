str_thisfunction = 'test_dolphin_2_collect_0';

disp(sprintf(' %% running %s',str_thisfunction));
%%%%%%%%;
setup_local;
nf=0;

n_x = 2;
dt_avg = 0.257;
n_j_factor = 1.02;
ignore_factor = 0.02;
T_ini = 0.00;
T_max = 10.5;
X_log_amplitude = 0.00;
flag_regularize_eccentricity_simultaneous = 1;

dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PAD');
dir_dolphin_mat = sprintf('%s/dir_dolphin_mat',dir_trunk);
%%%%;
n_i_ = [ 35 ; 70 ; 140 ]; n_n_i = numel(n_i_);
B_log_amplitude_ = [ -6 ; 0 ; +6 ]; n_B_log_amplitude = numel(B_log_amplitude_);
C_log_amplitude_ = [ -6 ; 0 ; +6 ]; n_C_log_amplitude = numel(C_log_amplitude_);
MaxFunEvals_use_simultaneous_ = [0 ; 16 ; 64]; n_MaxFunEvals_use_simultaneous = numel(MaxFunEvals_use_simultaneous_);
rseed_ = [1:32]; n_rseed = numel(rseed_);
dim_ = flip([n_n_i,n_B_log_amplitude,n_C_log_amplitude,n_MaxFunEvals_use_simultaneous,n_rseed]);
n_l = prod(dim_);
%%%%;
[rseed_l_,MaxFunEvals_use_simultaneous_l_,C_log_amplitude_l_,B_log_amplitude_l_,n_i_l_] = ndgrid(rseed_,MaxFunEvals_use_simultaneous_,C_log_amplitude_,B_log_amplitude_,n_i_);
rseed_l_ = rseed_l_(:);
MaxFunEvals_use_simultaneous_l_ = MaxFunEvals_use_simultaneous_l_(:);
C_log_amplitude_l_ = C_log_amplitude_l_(:);
B_log_amplitude_l_ = B_log_amplitude_l_(:);
n_i_l_ = n_i_l_(:);
%%%%;

n_t_sum_l_ = zeros(n_l,1); n_j_sum_l_ = zeros(n_l,1);
snr_A_vs_BC_l_ = zeros(n_l,1);
rfnorm_A_l_ = zeros(n_l,1); corr_A_l_ = zeros(n_l,1);
rfnorm_B_l_ = zeros(n_l,1); corr_B_l_ = zeros(n_l,1);
rfnorm_B_inv_l_ = zeros(n_l,1); corr_B_inv_l_ = zeros(n_l,1);
rfnorm_C_l_ = zeros(n_l,1); corr_C_l_ = zeros(n_l,1);
rfnorm_C_inv_l_ = zeros(n_l,1); corr_C_inv_l_ = zeros(n_l,1);
A_tru_eig_xl__ = zeros(n_x,n_l);
nl=0;
for nn_i=0:n_n_i-1;
n_i = n_i_(1+nn_i);
for nB_log_amplitude=0:n_B_log_amplitude-1;
B_log_amplitude = B_log_amplitude_(1+nB_log_amplitude);
for nC_log_amplitude=0:n_C_log_amplitude-1;
C_log_amplitude = C_log_amplitude_(1+nC_log_amplitude);
for nMaxFunEvals_use_simultaneous=0:n_MaxFunEvals_use_simultaneous-1;
MaxFunEvals_use_simultaneous = MaxFunEvals_use_simultaneous_(1+nMaxFunEvals_use_simultaneous);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
%%%%;
if (mod(nl,128)==0); disp(sprintf(' %% nl %.5d/%.5d',nl,n_l)); end;
str_infix ...
= ...
test_dolphin_infix_0( ...
 dt_avg ...
,n_i ...
,n_j_factor ...
,ignore_factor ...
,T_ini ...
,T_max ...
,B_log_amplitude ...
,C_log_amplitude ...
,X_log_amplitude ...
,MaxFunEvals_use_simultaneous ...
,flag_regularize_eccentricity_simultaneous ...
,rseed ...
);
fname_mat = sprintf('%s/test_dolphin_%s.mat',dir_dolphin_mat,str_infix);
if ~exist(fname_mat,'file'); disp(sprintf(' %% %s not found',fname_mat)); end;
if  exist(fname_mat,'file');
tmp_mat_ = load(fname_mat);
n_t_sum_l_(1+nl) = sum(tmp_mat_.n_t_i_); n_j_sum_l_(1+nl) = sum(tmp_mat_.n_j_i_);
snr_A_vs_BC_l_(1+nl) = sqrt(tmp_mat_.APdt_L2_avg)./max(1e-12,sqrt(tmp_mat_.BidW_L2_avg + tmp_mat_.Cinv_L2_avg));
snr_AB_vs_C_l_(1+nl) = sqrt(tmp_mat_.APdt_L2_avg + tmp_mat_.BidW_L2_avg)./max(1e-12,sqrt(tmp_mat_.Cinv_L2_avg));
A_tru_eig_xl__(:,1+nl) = eig(tmp_mat_.A_tru_xx__);
rfnorm_A_l_(1+nl) = tmp_mat_.rfnorm_A; corr_A_l_(1+nl) = tmp_mat_.corr_A;
rfnorm_B_l_(1+nl) = tmp_mat_.rfnorm_B; corr_B_l_(1+nl) = tmp_mat_.corr_B;
rfnorm_B_inv_l_(1+nl) = tmp_mat_.rfnorm_B_inv; corr_B_inv_l_(1+nl) = tmp_mat_.corr_B_inv;
rfnorm_C_l_(1+nl) = tmp_mat_.rfnorm_C; corr_C_l_(1+nl) = tmp_mat_.corr_C;
rfnorm_C_inv_l_(1+nl) = tmp_mat_.rfnorm_C_inv; corr_C_inv_l_(1+nl) = tmp_mat_.corr_C_inv;
nl=nl+1;
clear tmp_mat_;
end;%if  exist(fname_mat,'file');
%%%%;
end;%for nrseed=0:n_rseed-1;
end;%for nMaxFunEvals_use_simultaneous=0:n_MaxFunEvals_use_simultaneous-1;
end;%for nC_log_amplitude=0:n_C_log_amplitude-1;
end;%for nB_log_amplitude=0:n_B_log_amplitude-1;
end;%for nn_i=0:n_n_i-1;
%%%%%%%%;

figure(1+nf);nf=nf+1;clf;figbig;
p_row=2;p_col=4;np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_AB_vs_C_l_),log10(rfnorm_A_l_),16,C_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_AB_vs_C_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
xlim([-3,+3]);ylim([-3,+3]);grid on;
title('C_log_amplitude_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_AB_vs_C_l_),log10(rfnorm_A_l_),16,B_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_AB_vs_C_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
xlim([-3,+3]);ylim([-3,+3]);grid on;
title('B_log_amplitude_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),log10(rfnorm_A_l_),16,C_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
xlim([-3,+3]);ylim([-3,+3]);grid on;
title('C_log_amplitude_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),log10(rfnorm_A_l_),16,B_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
xlim([-3,+3]);ylim([-3,+3]);grid on;
title('B_log_amplitude_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_AB_vs_C_l_),corr_A_l_,16,C_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_AB_vs_C_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
xlim([-3,+3]);ylim([-1,+1]);grid on;
title('C_log_amplitude_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_AB_vs_C_l_),corr_A_l_,16,B_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_AB_vs_C_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
xlim([-3,+3]);ylim([-1,+1]);grid on;
title('B_log_amplitude_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),corr_A_l_,16,C_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
xlim([-3,+3]);ylim([-1,+1]);grid on;
title('C_log_amplitude_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),corr_A_l_,16,B_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
xlim([-3,+3]);ylim([-1,+1]);grid on;
title('B_log_amplitude_l_','Interpreter','none');
%%%%;


figure(1+nf);nf=nf+1;clf;figbig;
p_row=2;p_col=4;np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_AB_vs_C_l_),log10(rfnorm_A_l_),16,MaxFunEvals_use_simultaneous_l_,'filled');fig80s;
xlabel('log10(snr_AB_vs_C_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
xlim([-3,+3]);ylim([-3,+3]);grid on;
title('MaxFunEvals_use_simultaneous_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_AB_vs_C_l_),log10(rfnorm_A_l_),16,n_i_l_,'filled');fig80s;
xlabel('log10(snr_AB_vs_C_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
xlim([-3,+3]);ylim([-3,+3]);grid on;
title('n_i_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),log10(rfnorm_A_l_),16,MaxFunEvals_use_simultaneous_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
xlim([-3,+3]);ylim([-3,+3]);grid on;
title('MaxFunEvals_use_simultaneous_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),log10(rfnorm_A_l_),16,n_i_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
xlim([-3,+3]);ylim([-3,+3]);grid on;
title('n_i_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_AB_vs_C_l_),corr_A_l_,16,MaxFunEvals_use_simultaneous_l_,'filled');fig80s;
xlabel('log10(snr_AB_vs_C_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
xlim([-3,+3]);ylim([-1,+1]);grid on;
title('MaxFunEvals_use_simultaneous_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_AB_vs_C_l_),corr_A_l_,16,n_i_l_,'filled');fig80s;
xlabel('log10(snr_AB_vs_C_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
xlim([-3,+3]);ylim([-1,+1]);grid on;
title('n_i_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),corr_A_l_,16,MaxFunEvals_use_simultaneous_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
xlim([-3,+3]);ylim([-1,+1]);grid on;
title('MaxFunEvals_use_simultaneous_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),corr_A_l_,16,n_i_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
xlim([-3,+3]);ylim([-1,+1]);grid on;
title('n_i_l_','Interpreter','none');
%%%%;





