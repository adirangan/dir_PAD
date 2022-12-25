str_thisfunction = 'test_dolphin_3d_4_collect_0';

disp(sprintf(' %% running %s',str_thisfunction));
%%%%%%%%;
setup_local;
flag_replot = 1; 
nf=0;

n_x = 3;
dt_avg = 0.257;
n_j_factor = 1.02;
ignore_factor = 0.02;
T_ini = 0.00;
T_max = 10.5;
X_log_amplitude = 0.00;
flag_regularize_eccentricity_simultaneous = 1;

dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PAD');
dir_dolphin_mat = sprintf('%s/dir_dolphin_3d_mat',dir_trunk);
%%%%;
setup_local;
dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PAD');
n_i = 140;
n_iteration_BtBn = 32;
MaxFunEvals_use_BtBn = 1024*8;
MaxFunEvals_use_simultaneous = 0;
B_log_amplitude_ = [ -6:+3:+6 ]; n_B_log_amplitude = numel(B_log_amplitude_);
C_log_amplitude_ = [ -6:+3:+6 ]; n_C_log_amplitude = numel(C_log_amplitude_);
rseed_max = 1024;
rseed_ = [1:rseed_max]; n_rseed = numel(rseed_);

fname_collect_mat = sprintf('%s/test_dolphin_3d_4_collect_0.mat',dir_dolphin_mat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~exist(fname_collect_mat,'file');
disp(sprintf(' %% %s not found, creating',fname_collect_mat));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
dim_ = flip([n_B_log_amplitude,n_C_log_amplitude,n_rseed]);
n_l = prod(dim_);
%%%%;
[rseed_l_,C_log_amplitude_l_,B_log_amplitude_l_] = ndgrid(rseed_,C_log_amplitude_,B_log_amplitude_);
rseed_l_ = rseed_l_(:);
C_log_amplitude_l_ = C_log_amplitude_l_(:);
B_log_amplitude_l_ = B_log_amplitude_l_(:);
%%%%;
n_t_sum_l_ = zeros(n_l,1); n_j_sum_l_ = zeros(n_l,1);
snr_A_vs_BC_l_ = zeros(n_l,1);
snr_AB_vs_C_l_ = zeros(n_l,1);
snr_A_vs_B_l_ = zeros(n_l,1);
snr_A_vs_C_l_ = zeros(n_l,1);
fnorm_A_l_ = zeros(n_l,1);
fnorm_B_l_ = zeros(n_l,1);
fnorm_C_l_ = zeros(n_l,1);
rfnorm_A_l_ = zeros(n_l,1); corr_A_l_ = zeros(n_l,1);
rfnorm_B_l_ = zeros(n_l,1); corr_B_l_ = zeros(n_l,1);
rfnorm_B_inv_l_ = zeros(n_l,1); corr_B_inv_l_ = zeros(n_l,1);
rfnorm_C_l_ = zeros(n_l,1); corr_C_l_ = zeros(n_l,1);
rfnorm_C_inv_l_ = zeros(n_l,1); corr_C_inv_l_ = zeros(n_l,1);
A_tru_eig_xl__ = zeros(n_x,n_l);
nl=0;
for nB_log_amplitude=0:n_B_log_amplitude-1;
B_log_amplitude = B_log_amplitude_(1+nB_log_amplitude);
for nC_log_amplitude=0:n_C_log_amplitude-1;
C_log_amplitude = C_log_amplitude_(1+nC_log_amplitude);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
%%%%;
if (mod(nl,1024)==0); disp(sprintf(' %% nl %.5d/%.5d',nl,n_l)); end;
dir_infix ...
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
,n_iteration_BtBn ...
,MaxFunEvals_use_BtBn ...
,MaxFunEvals_use_simultaneous ...
,flag_regularize_eccentricity_simultaneous ...
,[] ...
);
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
,n_iteration_BtBn ...
,MaxFunEvals_use_BtBn ...
,MaxFunEvals_use_simultaneous ...
,flag_regularize_eccentricity_simultaneous ...
,rseed ...
);
fname_mat = sprintf('%s/dir_test_dolphin_3d_%s/test_dolphin_3d_%s.mat',dir_dolphin_mat,dir_infix,str_infix);
if ~exist(fname_mat,'file'); disp(sprintf(' %% %s not found',fname_mat)); end;
if  exist(fname_mat,'file');
tmp_mat_ = load(fname_mat);
n_t_sum_l_(1+nl) = sum(tmp_mat_.n_t_i_); n_j_sum_l_(1+nl) = sum(tmp_mat_.n_j_i_);
snr_A_vs_BC_l_(1+nl) = sqrt(tmp_mat_.APdt_L2_avg)./max(1e-12,sqrt(tmp_mat_.BidW_L2_avg + tmp_mat_.Cinv_L2_avg));
snr_AB_vs_C_l_(1+nl) = sqrt(tmp_mat_.APdt_L2_avg + tmp_mat_.BidW_L2_avg)./max(1e-12,sqrt(tmp_mat_.Cinv_L2_avg));
snr_A_vs_B_l_(1+nl) = sqrt(tmp_mat_.APdt_L2_avg)./max(1e-12,sqrt(tmp_mat_.BidW_L2_avg));
snr_A_vs_C_l_(1+nl) = sqrt(tmp_mat_.APdt_L2_avg)./max(1e-12,sqrt(tmp_mat_.Cinv_L2_avg));
A_tru_eig_xl__(:,1+nl) = eig(tmp_mat_.A_tru_xx__);
fnorm_A_l_(1+nl) = fnorm(tmp_mat_.A_tru_xx__); fnorm_B_l_(1+nl) = fnorm(sqrtm(tmp_mat_.BtBn_tru_inv_xx__)); fnorm_C_l_(1+nl) = fnorm(sqrtm(tmp_mat_.CtCn_tru_inv_xx__));
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
end;%for nC_log_amplitude=0:n_C_log_amplitude-1;
end;%for nB_log_amplitude=0:n_B_log_amplitude-1;
%%%%%%%%;
save(fname_collect_mat ...
     ,'dim_' ...
     ,'n_l' ...
     ,'rseed_l_' ...
     ,'C_log_amplitude_l_' ...
     ,'B_log_amplitude_l_' ...
     ,'n_t_sum_l_' ...
     ,'n_j_sum_l_' ...
     ,'snr_A_vs_BC_l_' ...
     ,'snr_AB_vs_C_l_' ...
     ,'snr_A_vs_B_l_' ...
     ,'snr_A_vs_C_l_' ...
     ,'fnorm_A_l_' ...
     ,'fnorm_B_l_' ...
     ,'fnorm_C_l_' ...
     ,'rfnorm_A_l_' ...
     ,'rfnorm_B_l_' ...
     ,'rfnorm_B_inv_l_' ...
     ,'rfnorm_C_l_' ...
     ,'rfnorm_C_inv_l_' ...
     ,'corr_A_l_' ...
     ,'corr_B_l_' ...
     ,'corr_B_inv_l_' ...
     ,'corr_C_l_' ...
     ,'corr_C_inv_l_' ...
     ,'A_tru_eig_xl__' ...
     );
%%%%%%%%;
end;%if ~exist(fname_collect_mat,'file');
if  exist(fname_collect_mat,'file');
disp(sprintf(' %% %s found, not creating',fname_collect_mat));
load(fname_collect_mat);
end;%if  exist(fname_collect_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% examine relationship between snr_A_vs_BC and snr_AB_vs_C and error. ;
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
xlim([-3,+1]);ylim([-3,+3]);grid on;
title('C_log_amplitude_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),log10(rfnorm_A_l_),16,B_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
xlim([-3,+1]);ylim([-3,+3]);grid on;
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
xlim([-3,+1]);ylim([-1,+1]);grid on;
title('C_log_amplitude_l_','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),corr_A_l_,16,B_log_amplitude_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
xlim([-3,+1]);ylim([-1,+1]);grid on;
title('B_log_amplitude_l_','Interpreter','none');
%%%%;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% plot error (corr_A) as function of fnorm_C for various fnorm_B. ;
%%%%%%%%;
l10_fnorm_C_lim_ = [-1.5,+1.5]; n_l10_fnorm_C_bin = 1+16; 
l10_fnorm_B_lim_ = [-2.0,+2.0]; n_l10_fnorm_B_bin = 1+4;
l10_rfnorm_A_lim_ = [-1.5,+1.5]; n_l10_rfnorm_A_bin = 1+16;
corr_A_lim_ = [-1.0,+1.0]; n_corr_A_bin = 1+20;
%%%%;
l10_fnorm_C_tick_ = 1:2:n_l10_fnorm_C_bin; l10_fnorm_C_bin_ = linspace(min(l10_fnorm_C_lim_),max(l10_fnorm_C_lim_),n_l10_fnorm_C_bin);
l10_fnorm_B_tick_ = 1:2:n_l10_fnorm_B_bin; l10_fnorm_B_bin_ = linspace(min(l10_fnorm_B_lim_),max(l10_fnorm_B_lim_),n_l10_fnorm_B_bin);
l10_rfnorm_A_tick_ = 1:2:n_l10_rfnorm_A_bin; l10_rfnorm_A_bin_ = linspace(min(l10_rfnorm_A_lim_),max(l10_rfnorm_A_lim_),n_l10_rfnorm_A_bin);
corr_A_tick_ = 1:2:n_corr_A_bin; corr_A_bin_ = linspace(min(corr_A_lim_),max(corr_A_lim_),n_corr_A_bin);
%%%%;
[h3d_CBA___] = hist3d_0(log10(fnorm_C_l_),log10(fnorm_B_l_),corr_A_l_,n_l10_fnorm_C_bin,n_l10_fnorm_B_bin,n_corr_A_bin,l10_fnorm_C_lim_,l10_fnorm_B_lim_,corr_A_lim_);
h3d_CBA___ = bsxfun(@rdivide,h3d_CBA___,max(1,sum(h3d_CBA___,3)));
%%%%;
figure(1+nf);nf=nf+1;clf;figbig; fig80s;
%p_row = 2; p_col = ceil(n_l10_fnorm_B_bin/p_row); np=0;
p_row = 2; p_col = ceil(2/p_row); np=0;
%for np=0:n_l10_fnorm_B_bin-1;
for np=0:p_row*p_col-1;
nl10_fnorm_B_bin = floor(n_l10_fnorm_B_bin/2)+np;
subplot(p_row,p_col,1+np); cla;
h2d_CA__ = squeeze(h3d_CBA___(:,1+nl10_fnorm_B_bin,:));
imagesc(transpose(h2d_CA__));
set(gca,'ydir','normal'); %axis image;
set(gca,'XTick',l10_fnorm_C_tick_,'XTickLabel',l10_fnorm_C_bin_(l10_fnorm_C_tick_));
set(gca,'YTick',corr_A_tick_,'YTickLabel',corr_A_bin_(corr_A_tick_));
xlabel('l10_fnorm_C_lim_','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
title(sprintf('l10_fnorm_B_lim_ %0.2f',min(l10_fnorm_B_lim_) + nl10_fnorm_B_bin*diff(l10_fnorm_B_lim_)/max(1,n_l10_fnorm_B_bin-1)),'Interpreter','none');
end;%for np=0:n_l10_fnorm_B_bin-1;
%%%%;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% examing relationship between fnorm_C, fnorm_B and error. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
p_row=1;p_col=1;np=0;
l10_fnorm_C_lim_ = [-2.0,+2.0];
l10_fnorm_B_lim_ = [-2.0,+2.0];
l10_rfnorm_A_lim_ = [-1.5,+1.5];
corr_A_lim_ = [-1.1,+1.1];
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter3(log10(fnorm_C_l_),log10(fnorm_B_l_),log10(rfnorm_A_l_),16,corr_A_l_,'filled');fig80s;
xlabel('log10(fnorm_C_l_)','Interpreter','none');
ylabel('log10(fnorm_B_l_)','Interpreter','none');
zlabel('log10(rfnorm_A_l_)','Interpreter','none');
title('corr_A_l_','Interpreter','none');
xlim(l10_fnorm_C_lim_); ylim(l10_fnorm_B_lim_); zlim(l10_rfnorm_A_lim_); grid on; axis vis3d;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
p_row=1;p_col=1;np=0;
l10_fnorm_C_lim_ = [-2.0,+2.0];
l10_fnorm_B_lim_ = [-2.0,+2.0];
l10_rfnorm_A_lim_ = [-1.5,+1.5];
corr_A_lim_ = [-1.1,+1.1];
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter3(log10(fnorm_C_l_),log10(fnorm_B_l_),corr_A_l_,16,log10(rfnorm_A_l_),'filled');fig80s;
xlabel('log10(fnorm_C_l_)','Interpreter','none');
ylabel('log10(fnorm_B_l_)','Interpreter','none');
zlabel('corr_A_l_','Interpreter','none');
title('log10(rfnorm_A_l_)','Interpreter','none');
xlim(l10_fnorm_C_lim_); ylim(l10_fnorm_B_lim_); zlim(corr_A_lim_); grid on; axis vis3d;
%%%%%%%%;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% examing relationship between snr_A_vs_C, snr_A_vs_B and error. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
p_row=1;p_col=1;np=0;
l10_snr_A_vs_C_lim_ = [-3.5,+3.5];
l10_snr_A_vs_B_lim_ = [-1.5,+1.5];
l10_rfnorm_A_lim_ = [-1.5,+1.5];
corr_A_lim_ = [-1.1,+1.1];
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter3(log10(snr_A_vs_C_l_),log10(snr_A_vs_B_l_),log10(rfnorm_A_l_),16,corr_A_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_C_l_)','Interpreter','none');
ylabel('log10(snr_A_vs_B_l_)','Interpreter','none');
zlabel('log10(rfnorm_A_l_)','Interpreter','none');
title('corr_A_l_','Interpreter','none');
xlim(l10_snr_A_vs_C_lim_); ylim(l10_snr_A_vs_B_lim_); zlim(l10_rfnorm_A_lim_); grid on; axis vis3d;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
p_row=1;p_col=1;np=0;
l10_snr_A_vs_C_lim_ = [-3.5,+3.5];
l10_snr_A_vs_B_lim_ = [-1.5,+1.5];
l10_rfnorm_A_lim_ = [-1.5,+1.5];
corr_A_lim_ = [-1.1,+1.1];
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter3(log10(snr_A_vs_C_l_),log10(snr_A_vs_B_l_),corr_A_l_,16,log10(rfnorm_A_l_),'filled');fig80s;
xlabel('log10(snr_A_vs_C_l_)','Interpreter','none');
ylabel('log10(snr_A_vs_B_l_)','Interpreter','none');
zlabel('corr_A_l_','Interpreter','none');
title('log10(rfnorm_A_l_)','Interpreter','none');
xlim(l10_snr_A_vs_C_lim_); ylim(l10_snr_A_vs_B_lim_); zlim(corr_A_lim_); grid on; axis vis3d;
%%%%%%%%;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% examine relationship between snr_A_vs_BC and snr_AB_vs_C and error. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row=2;p_col=3;np=0;
l10_snr_A_vs_BC_lim_ = [-3.5,+1];
l10_snr_A_vs_B_lim_ = [-1.5,+1.5];
l10_rfnorm_A_lim_ = [-1.5,+1.5];
corr_A_lim_ = [-1.1,+1.1];
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),log10(rfnorm_A_l_),16,log10(snr_A_vs_B_l_),'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
title('log10(snr_A_vs_B_l_)','Interpreter','none');
xlim(l10_snr_A_vs_BC_lim_); ylim(l10_rfnorm_A_lim_); grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_B_l_),log10(rfnorm_A_l_),16,log10(snr_A_vs_BC_l_),'filled');fig80s;
xlabel('log10(snr_A_vs_B_l_)','Interpreter','none');
ylabel('log10(rfnorm_A_l_)','Interpreter','none');
title('log10(snr_A_vs_BC_l_)','Interpreter','none');
xlim(l10_snr_A_vs_B_lim_); ylim(l10_rfnorm_A_lim_); grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),log10(snr_A_vs_B_l_),16,log10(rfnorm_A_l_),'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('log10(snr_A_vs_B_l_)','Interpreter','none');
title('log10(rfnorm_A_l_)','Interpreter','none');
xlim(l10_snr_A_vs_BC_lim_); ylim(l10_snr_A_vs_B_lim_); grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),corr_A_l_,16,log10(snr_A_vs_B_l_),'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
title('log10(snr_A_vs_B_l_)','Interpreter','none');
xlim(l10_snr_A_vs_BC_lim_); ylim(corr_A_lim_); grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_B_l_),corr_A_l_,16,log10(snr_A_vs_BC_l_),'filled');fig80s;
xlabel('log10(snr_A_vs_B_l_)','Interpreter','none');
ylabel('corr_A_l_','Interpreter','none');
title('log10(snr_A_vs_BC_l_)','Interpreter','none');
xlim(l10_snr_A_vs_B_lim_); ylim(corr_A_lim_); grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
scatter(log10(snr_A_vs_BC_l_),log10(snr_A_vs_B_l_),16,corr_A_l_,'filled');fig80s;
xlabel('log10(snr_A_vs_BC_l_)','Interpreter','none');
ylabel('log10(snr_A_vs_B_l_)','Interpreter','none');
title('corr_A_l_','Interpreter','none');
xlim(l10_snr_A_vs_BC_lim_); ylim(l10_snr_A_vs_B_lim_); grid on;
%%%%%%%%;
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% plot histogram of corr_A for each fnorm_C and a few fnorm_B. ; 
%%%%%%%%;
l10_fnorm_B_lo_ = [-0.30,+0.30];
index_fnorm_B_lo_ = efind( (log10(fnorm_B_l_)>=min(l10_fnorm_B_lo_)) & (log10(fnorm_B_l_)<=max(l10_fnorm_B_lo_)) );
l10_fnorm_B_me_ = [+0.30,+0.90];
index_fnorm_B_me_ = efind( (log10(fnorm_B_l_)>=min(l10_fnorm_B_me_)) & (log10(fnorm_B_l_)<=max(l10_fnorm_B_me_)) );
l10_fnorm_B_hi_ = [+0.50,+2.10];
index_fnorm_B_hi_ = efind( (log10(fnorm_B_l_)>=min(l10_fnorm_B_hi_)) & (log10(fnorm_B_l_)<=max(l10_fnorm_B_hi_)) );
l10_fnorm_C_lim_ = [-1.5,+1.5]; n_l10_fnorm_C_bin = 1+7; 
l10_fnorm_C_bin_ = linspace(min(l10_fnorm_C_lim_),max(l10_fnorm_C_lim_),n_l10_fnorm_C_bin);
l10_fnorm_C_tick_ = 1:1:n_l10_fnorm_C_bin; 
l10_fnorm_C_ticklabel_ = num2str(transpose(l10_fnorm_C_bin_(l10_fnorm_C_tick_)),'%+.1f');
corr_A_lim_ = [-1.00,+1.00]; n_corr_A_bin = 1+40;
corr_A_bin_ = linspace(min(corr_A_lim_),max(corr_A_lim_),n_corr_A_bin);
ij_corr_A_bin_use_ = ceil(n_corr_A_bin/2):n_corr_A_bin;
corr_A_bin_use_ = corr_A_bin_(ij_corr_A_bin_use_); n_corr_A_bin_use = numel(corr_A_bin_use_);
corr_A_bin_use_lim_ = [min(corr_A_bin_use_),max(corr_A_bin_use_)];
corr_A_tick_ = 1:2:n_corr_A_bin_use;
h2d_lo_AC__ = hist2d_0(log10(fnorm_C_l_(1+index_fnorm_B_lo_)),corr_A_l_(1+index_fnorm_B_lo_),n_l10_fnorm_C_bin,n_corr_A_bin,l10_fnorm_C_lim_,corr_A_lim_);
h2d_lo_AC__ = bsxfun(@rdivide,h2d_lo_AC__,sum(h2d_lo_AC__,1)) / max(1e-12,mean(diff(corr_A_bin_)));
h2d_me_AC__ = hist2d_0(log10(fnorm_C_l_(1+index_fnorm_B_me_)),corr_A_l_(1+index_fnorm_B_me_),n_l10_fnorm_C_bin,n_corr_A_bin,l10_fnorm_C_lim_,corr_A_lim_);
h2d_me_AC__ = bsxfun(@rdivide,h2d_me_AC__,sum(h2d_me_AC__,1)) / max(1e-12,mean(diff(corr_A_bin_)));
h2d_hi_AC__ = hist2d_0(log10(fnorm_C_l_(1+index_fnorm_B_hi_)),corr_A_l_(1+index_fnorm_B_hi_),n_l10_fnorm_C_bin,n_corr_A_bin,l10_fnorm_C_lim_,corr_A_lim_);
h2d_hi_AC__ = bsxfun(@rdivide,h2d_hi_AC__,sum(h2d_hi_AC__,1)) / max(1e-12,mean(diff(corr_A_bin_)));
%%%%;
corr_A_avg_lo_avg_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_me_avg_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_hi_avg_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_lo_p15_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_me_p15_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_hi_p15_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_lo_p50_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_me_p50_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_hi_p50_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_lo_p85_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_me_p85_b_ = zeros(n_l10_fnorm_C_bin,1);
corr_A_avg_hi_p85_b_ = zeros(n_l10_fnorm_C_bin,1);
for nl10_fnorm_C_bin=0:n_l10_fnorm_C_bin-1;
edge_0 = -Inf; if (nl10_fnorm_C_bin> 0); edge_0 = 0.5*(l10_fnorm_C_bin_(1+nl10_fnorm_C_bin-1) + l10_fnorm_C_bin_(1+nl10_fnorm_C_bin-0)); end;
edge_1 = +Inf; if (nl10_fnorm_C_bin< n_l10_fnorm_C_bin-1); edge_1 = 0.5*(l10_fnorm_C_bin_(1+nl10_fnorm_C_bin+1) + l10_fnorm_C_bin_(1+nl10_fnorm_C_bin+0)); end;
tmp_index_lo_ = efind( (log10(fnorm_C_l_(1+index_fnorm_B_lo_))>=edge_0) & (log10(fnorm_C_l_(1+index_fnorm_B_lo_))<=edge_1) );
corr_A_avg_lo_avg_b_(1+nl10_fnorm_C_bin) = mean(corr_A_l_(1+index_fnorm_B_lo_(1+tmp_index_lo_)),'omitnan');
corr_A_avg_lo_p15_b_(1+nl10_fnorm_C_bin) = prctile(corr_A_l_(1+index_fnorm_B_lo_(1+tmp_index_lo_)),15);
corr_A_avg_lo_p50_b_(1+nl10_fnorm_C_bin) = prctile(corr_A_l_(1+index_fnorm_B_lo_(1+tmp_index_lo_)),50);
corr_A_avg_lo_p85_b_(1+nl10_fnorm_C_bin) = prctile(corr_A_l_(1+index_fnorm_B_lo_(1+tmp_index_lo_)),85);
tmp_index_me_ = efind( (log10(fnorm_C_l_(1+index_fnorm_B_me_))>=edge_0) & (log10(fnorm_C_l_(1+index_fnorm_B_me_))<=edge_1) );
corr_A_avg_me_avg_b_(1+nl10_fnorm_C_bin) = mean(corr_A_l_(1+index_fnorm_B_me_(1+tmp_index_me_)),'omitnan');
corr_A_avg_me_p15_b_(1+nl10_fnorm_C_bin) = prctile(corr_A_l_(1+index_fnorm_B_me_(1+tmp_index_me_)),15);
corr_A_avg_me_p50_b_(1+nl10_fnorm_C_bin) = prctile(corr_A_l_(1+index_fnorm_B_me_(1+tmp_index_me_)),50);
corr_A_avg_me_p85_b_(1+nl10_fnorm_C_bin) = prctile(corr_A_l_(1+index_fnorm_B_me_(1+tmp_index_me_)),85);
tmp_index_hi_ = efind( (log10(fnorm_C_l_(1+index_fnorm_B_hi_))>=edge_0) & (log10(fnorm_C_l_(1+index_fnorm_B_hi_))<=edge_1) );
corr_A_avg_hi_avg_b_(1+nl10_fnorm_C_bin) = mean(corr_A_l_(1+index_fnorm_B_hi_(1+tmp_index_hi_)),'omitnan');
corr_A_avg_hi_p15_b_(1+nl10_fnorm_C_bin) = prctile(corr_A_l_(1+index_fnorm_B_hi_(1+tmp_index_hi_)),15);
corr_A_avg_hi_p50_b_(1+nl10_fnorm_C_bin) = prctile(corr_A_l_(1+index_fnorm_B_hi_(1+tmp_index_hi_)),50);
corr_A_avg_hi_p85_b_(1+nl10_fnorm_C_bin) = prctile(corr_A_l_(1+index_fnorm_B_hi_(1+tmp_index_hi_)),85);
disp(sprintf(' %% lo %d me %d hi %d',numel(tmp_index_lo_),numel(tmp_index_me_),numel(tmp_index_hi_)));
end;%for nl10_fnorm_C_bin=0:n_l10_fnorm_C_bin-1;
%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,2048*0.75,512*0.75]); colormap('hot');
p_row = 1; p_col = 3; np=0; lclim_ = [-2,+2];
fontsize_use = 12; linewidth_big = 4; linewidth_sml = 2; markersize_use = 16;
%%%%;
subplot(p_row,p_col,1+np);np=np+1; cla;
hold on;
imagesc(log2(0+h2d_lo_AC__(ij_corr_A_bin_use_,:)),lclim_); tmp_c_ = colorbar; set(tmp_c_,'TickLength',0);
plot(1.0 + [0:n_l10_fnorm_C_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_lo_p15_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_sml,'MarkerSize',markersize_use);
plot(1.0 + [0:n_l10_fnorm_C_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_lo_p50_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_big,'MarkerSize',markersize_use);
plot(1.0 + [0:n_l10_fnorm_C_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_lo_p85_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_sml,'MarkerSize',markersize_use);
hold off;
xlim([0.5,0.5+n_l10_fnorm_C_bin]); ylim([0.5,0.5+n_corr_A_bin_use]);
set(gca,'ydir','normal'); %axis image;
set(gca,'XTick',l10_fnorm_C_tick_,'XTickLabel',l10_fnorm_C_ticklabel_); xtickangle(90);
set(gca,'YTick',corr_A_tick_,'YTickLabel',corr_A_bin_use_(corr_A_tick_));
xlabel('log10(magnitude of C)','Interpreter','none');
ylabel('correlation with A','Interpreter','none');
set(gca,'TickLength',[0,0]);
set(gca,'FontSize',fontsize_use);
%title(sprintf('l10_fnorm_B_lim_ %0.2f,%0.2f',l10_fnorm_B_lo_),'Interpreter','none');
title(sprintf('B small'),'Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1; cla;
hold on;
imagesc(log2(0+h2d_me_AC__(ij_corr_A_bin_use_,:)),lclim_); tmp_c_ = colorbar; set(tmp_c_,'TickLength',0);
plot(1.0 + [0:n_l10_fnorm_C_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_me_p15_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_sml,'MarkerSize',markersize_use);
plot(1.0 + [0:n_l10_fnorm_C_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_me_p50_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_big,'MarkerSize',markersize_use);
plot(1.0 + [0:n_l10_fnorm_C_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_me_p85_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_sml,'MarkerSize',markersize_use);
hold off;
xlim([0.5,0.5+n_l10_fnorm_C_bin]); ylim([0.5,0.5+n_corr_A_bin_use]);
set(gca,'ydir','normal'); %axis image;
set(gca,'XTick',l10_fnorm_C_tick_,'XTickLabel',l10_fnorm_C_ticklabel_); xtickangle(90);
set(gca,'YTick',corr_A_tick_,'YTickLabel',corr_A_bin_use_(corr_A_tick_));
xlabel('log10(magnitude of C)','Interpreter','none');
ylabel('correlation with A','Interpreter','none');
set(gca,'TickLength',[0,0]);
set(gca,'FontSize',fontsize_use);
%title(sprintf('l10_fnorm_B_lim_ %0.2f,%0.2f',l10_fnorm_B_me_),'Interpreter','none');
title(sprintf('B medium'),'Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1; cla;
hold on;
imagesc(log2(0+h2d_hi_AC__(ij_corr_A_bin_use_,:)),lclim_); tmp_c_ = colorbar; set(tmp_c_,'TickLength',0);
plot(1.0 + [0:n_l10_fnorm_C_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_hi_p15_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_sml,'MarkerSize',markersize_use);
plot(1.0 + [0:n_l10_fnorm_C_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_hi_p50_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_big,'MarkerSize',markersize_use);
plot(1.0 + [0:n_l10_fnorm_C_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_hi_p85_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_sml,'MarkerSize',markersize_use);
hold off;
xlim([0.5,0.5+n_l10_fnorm_C_bin]); ylim([0.5,0.5+n_corr_A_bin_use]);
set(gca,'ydir','normal'); %axis image;
set(gca,'XTick',l10_fnorm_C_tick_,'XTickLabel',l10_fnorm_C_ticklabel_); xtickangle(90);
set(gca,'YTick',corr_A_tick_,'YTickLabel',corr_A_bin_use_(corr_A_tick_));
xlabel('log10(magnitude of C)','Interpreter','none');
ylabel('correlation with A','Interpreter','none');
set(gca,'TickLength',[0,0]);
set(gca,'FontSize',fontsize_use);
%title(sprintf('l10_fnorm_B_lim_ %0.2f,%0.2f',l10_fnorm_B_hi_),'Interpreter','none');
title(sprintf('B large'),'Interpreter','none');
%%%%%%%%;
dir_dolphin_jpg = sprintf('%s/dir_dolphin_3d_jpg',dir_trunk);
if ~exist(dir_dolphin_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_dolphin_jpg)); mkdir(dir_dolphin_jpg); end;
fname_fig_pre = sprintf('%s/test_dolphin_3d_4_collect_0_FIGA',dir_dolphin_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% plot histogram of corr_A for each snr_A_vs_BC. ;
%%%%%%%%;
l10_snr_A_vs_BC_lim_ = [-0.0,+2.5]; n_l10_snr_A_vs_BC_bin = 1+7; 
l10_snr_A_vs_BC_bin_ = linspace(min(l10_snr_A_vs_BC_lim_),max(l10_snr_A_vs_BC_lim_),n_l10_snr_A_vs_BC_bin);
l10_snr_A_vs_BC_tick_ = 1:1:n_l10_snr_A_vs_BC_bin; 
l10_snr_A_vs_BC_ticklabel_ = num2str(transpose(l10_snr_A_vs_BC_bin_(l10_snr_A_vs_BC_tick_)),'%+.1f');
corr_A_lim_ = [-1.00,+1.00]; n_corr_A_bin = 1+40;
corr_A_bin_ = linspace(min(corr_A_lim_),max(corr_A_lim_),n_corr_A_bin);
ij_corr_A_bin_use_ = ceil(n_corr_A_bin/2):n_corr_A_bin;
%ij_corr_A_bin_use_ = 1:n_corr_A_bin;
corr_A_bin_use_ = corr_A_bin_(ij_corr_A_bin_use_); n_corr_A_bin_use = numel(corr_A_bin_use_);
corr_A_bin_use_lim_ = [min(corr_A_bin_use_),max(corr_A_bin_use_)];
corr_A_tick_ = 1:2:n_corr_A_bin_use;
h2d_snr_A_vs_BC__ = hist2d_0(-log10(snr_A_vs_BC_l_),corr_A_l_,n_l10_snr_A_vs_BC_bin,n_corr_A_bin,l10_snr_A_vs_BC_lim_,corr_A_lim_);
h2d_snr_A_vs_BC__ = bsxfun(@rdivide,h2d_snr_A_vs_BC__,sum(h2d_snr_A_vs_BC__,1)) / max(1e-12,mean(diff(corr_A_bin_)));
%%%%;
corr_A_avg_avg_b_ = zeros(n_l10_snr_A_vs_BC_bin,1);
corr_A_avg_p15_b_ = zeros(n_l10_snr_A_vs_BC_bin,1);
corr_A_avg_p50_b_ = zeros(n_l10_snr_A_vs_BC_bin,1);
corr_A_avg_p85_b_ = zeros(n_l10_snr_A_vs_BC_bin,1);
for nl10_snr_A_vs_BC_bin=0:n_l10_snr_A_vs_BC_bin-1;
edge_0 = -Inf; if (nl10_snr_A_vs_BC_bin> 0); edge_0 = 0.5*(l10_snr_A_vs_BC_bin_(1+nl10_snr_A_vs_BC_bin-1) + l10_snr_A_vs_BC_bin_(1+nl10_snr_A_vs_BC_bin-0)); end;
edge_1 = +Inf; if (nl10_snr_A_vs_BC_bin< n_l10_snr_A_vs_BC_bin-1); edge_1 = 0.5*(l10_snr_A_vs_BC_bin_(1+nl10_snr_A_vs_BC_bin+1) + l10_snr_A_vs_BC_bin_(1+nl10_snr_A_vs_BC_bin+0)); end;
tmp_index_ = efind( (-log10(snr_A_vs_BC_l_)>=edge_0) & (-log10(snr_A_vs_BC_l_)<=edge_1) );
corr_A_avg_avg_b_(1+nl10_snr_A_vs_BC_bin) = mean(corr_A_l_(1+tmp_index_),'omitnan');
corr_A_avg_p15_b_(1+nl10_snr_A_vs_BC_bin) = prctile(corr_A_l_(1+tmp_index_),15);
corr_A_avg_p50_b_(1+nl10_snr_A_vs_BC_bin) = prctile(corr_A_l_(1+tmp_index_),50);
corr_A_avg_p85_b_(1+nl10_snr_A_vs_BC_bin) = prctile(corr_A_l_(1+tmp_index_),85);
end;%for nl10_snr_A_vs_BC_bin=0:n_l10_snr_A_vs_BC_bin-1;
%%%%;
figure(1+nf);nf=nf+1;clf;figsml; colormap('hot');
p_row = 1; p_col = 1; np=0; lclim_ = [-2,+2];
fontsize_use = 12; linewidth_big = 4; linewidth_sml = 2; markersize_use = 16;
%%%%;
subplot(p_row,p_col,1+np);np=np+1; cla;
hold on;
imagesc(log2(0+h2d_snr_A_vs_BC__(ij_corr_A_bin_use_,:)),lclim_); tmp_c_ = colorbar; set(tmp_c_,'TickLength',0);
plot(1.0 + [0:n_l10_snr_A_vs_BC_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_p15_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_sml,'MarkerSize',markersize_use);
plot(1.0 + [0:n_l10_snr_A_vs_BC_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_p50_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_big,'MarkerSize',markersize_use);
plot(1.0 + [0:n_l10_snr_A_vs_BC_bin-1],0.5 + n_corr_A_bin_use*(corr_A_avg_p85_b_-min(corr_A_bin_use_lim_))/diff(corr_A_bin_use_lim_),'c.-','LineWidth',linewidth_sml,'MarkerSize',markersize_use);
hold off;
xlim([0.5,0.5+n_l10_snr_A_vs_BC_bin]); ylim([0.5,0.5+n_corr_A_bin_use]);
set(gca,'ydir','normal'); %axis image;
set(gca,'XTick',l10_snr_A_vs_BC_tick_,'XTickLabel',l10_snr_A_vs_BC_ticklabel_); xtickangle(90);
set(gca,'YTick',corr_A_tick_,'YTickLabel',corr_A_bin_use_(corr_A_tick_));
xlabel('-log10 snr of A vs (B,C)','Interpreter','none');
ylabel('correlation with A','Interpreter','none');
set(gca,'TickLength',[0,0]);
set(gca,'FontSize',fontsize_use);
%title(sprintf('l10_fnorm_B_lim_ %0.2f,%0.2f',l10_fnorm_B_),'Interpreter','none');
title(sprintf(''),'Interpreter','none');
%%%%%%%%;
dir_dolphin_jpg = sprintf('%s/dir_dolphin_3d_jpg',dir_trunk);
if ~exist(dir_dolphin_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_dolphin_jpg)); mkdir(dir_dolphin_jpg); end;
fname_fig_pre = sprintf('%s/test_dolphin_3d_4_collect_0_FIGB',dir_dolphin_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
end;%if flag_disp;





