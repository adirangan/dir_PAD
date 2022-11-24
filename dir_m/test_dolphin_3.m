function ...
[ ...
 parameter_dolphin ...
] = ...
test_dolphin_3( ...
 parameter_dolphin ...
);
%%%%%%%%;
% Built from test_SDE_nlp_ijXaABYC_update_1.m. ;
% Using statistics from dolphin data. ;
% e.g.,: dt_avg = 0.257; %<-- matches d00 data (std 0.344) (var 0.118) ;
% n_i = 140; %<-- roughly matches d00 data. ;
% n_j_factor = 1.02; %<-- matches d00 data. ;
% ignore_factor = 0.02; <-- matches d00 data. ;
% T_ini =  0;
% T_max = 10.5;%T_max = 55*1;%T_max = 55*144;%T_max = 55; %<-- average T_max is roughly 10.12. ;
%%%%%%%%;
%{
  tmp_ = load('/data/rangan/dir_bcc/dir_dolphin/dir_mat/dolphin_aid00_dt_0.mat');
  u_aid_ = unique(tmp_.aid_); n_u_aid = numel(u_aid_);
  age_at__ = cell(n_u_aid,1);
  for nu_aid=0:n_u_aid-1;
  u_aid = u_aid_(1+nu_aid);
  tmp_index_ = efind(tmp_.aid_==u_aid);
  age_t_ = tmp_.age_(1+tmp_index_);
  age_at__{1+nu_aid} = age_t_;
  clear age_t_;
  end;%for nu_aid=0:n_u_aid-1;
  amax_a_ = zeros(n_u_aid,1);
  amin_a_ = zeros(n_u_aid,1);
  for nu_aid=0:n_u_aid-1;
  age_t_ = age_at__{1+nu_aid};
  amax_a_(1+nu_aid) = max(age_t_);
  amin_a_(1+nu_aid) = min(age_t_);
  clear age_t_;
  end;%for nu_aid=0:n_u_aid-1;
  adif_a_ = amax_a_ - amin_a_;
  disp(sprintf(' %% mean(adif_a_): %0.2f',mean(adif_a_)));
  % mean(adif_a_) is about 10. ;
  % in order to get about 5800 time-steps, ;
  % we need either: ;
  % 1 dolphin with T_max = 55*27, ;
  % 27 dolphins with T_max = 55, ; (Note that 55 is the max age of a single dolphin). ;
  % or 120 dolphins with T_max = 10. ; (Note that the observed number of dolphins was 136). ;
%}
%%%%%%%%;
str_thisfunction = 'test_dolphin_3';

if nargin<1;
disp(sprintf(' %% running %s',str_thisfunction));
%%%%%%%%;
setup_local;
dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PAD');
n_i = 140;
n_iteration_BtBn = 32;
MaxFunEvals_use_BtBn = 1024;
MaxFunEvals_use_simultaneous = 0;
B_log_amplitude_ = [ -6:+3:+6 ]; n_B_log_amplitude = numel(B_log_amplitude_);
C_log_amplitude_ = [ -6:+3:+6 ]; n_C_log_amplitude = numel(C_log_amplitude_);
rseed_ = [1:256]; n_rseed = numel(rseed_);
for nB_log_amplitude=0:n_B_log_amplitude-1;
B_log_amplitude = B_log_amplitude_(1+nB_log_amplitude);
for nC_log_amplitude=0:n_C_log_amplitude-1;
C_log_amplitude = C_log_amplitude_(1+nC_log_amplitude);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
%%%%;
parameter_dolphin = struct('type','parameter');
parameter_dolphin.dir_trunk = dir_trunk;
parameter_dolphin.n_i = n_i;
parameter_dolphin.B_log_amplitude = B_log_amplitude;
parameter_dolphin.C_log_amplitude = C_log_amplitude;
parameter_dolphin.n_iteration_BtBn = n_iteration_BtBn;
parameter_dolphin.MaxFunEvals_use_BtBn = MaxFunEvals_use_BtBn;
parameter_dolphin.MaxFunEvals_use_simultaneous = MaxFunEvals_use_simultaneous;
parameter_dolphin.rseed = rseed;
test_dolphin_3(parameter_dolphin);
%%%%;
end;%for nrseed=0:n_rseed-1;
end;%for nC_log_amplitude=0:n_C_log_amplitude-1;
end;%for nB_log_amplitude=0:n_B_log_amplitude-1;
%%%%%%%%;
disp(sprintf('returning'));return;
end;%if nargin<1;

%%%%%%%%;
na=0;
if (nargin<1+na); parameter_dolphin=[]; end; na=na+1;
%%%%%%%%;
if isempty(parameter_dolphin); parameter_dolphin = struct('type','parameter'); end;
if ~isfield(parameter_dolphin,'dir_trunk'); parameter_dolphin.dir_trunk = pwd; end;
if ~isfield(parameter_dolphin,'dir_dolphin_mat'); parameter_dolphin.dir_dolphin_mat = sprintf('%s/dir_dolphin_mat',parameter_dolphin.dir_trunk); end;
if ~isfield(parameter_dolphin,'tolerance_master'); parameter_dolphin.tolerance_master = 1e-6; end;
if ~isfield(parameter_dolphin,'flag_verbose'); parameter_dolphin.flag_verbose = 1; end;
if ~isfield(parameter_dolphin,'flag_disp'); parameter_dolphin.flag_disp = 0; end;
if ~isfield(parameter_dolphin,'dt_avg'); parameter_dolphin.dt_avg = 0.257; end;
if ~isfield(parameter_dolphin,'n_i'); parameter_dolphin.n_i = 140; end;
if ~isfield(parameter_dolphin,'n_j_factor'); parameter_dolphin.n_j_factor = 1.02; end;
if ~isfield(parameter_dolphin,'ignore_factor'); parameter_dolphin.ignore_factor = 0.02; end;
if ~isfield(parameter_dolphin,'T_ini'); parameter_dolphin.T_ini = 0.00; end;
if ~isfield(parameter_dolphin,'T_max'); parameter_dolphin.T_max = 10.5; end;
if ~isfield(parameter_dolphin,'B_log_amplitude'); parameter_dolphin.B_log_amplitude = 0.0; end;
if ~isfield(parameter_dolphin,'C_log_amplitude'); parameter_dolphin.C_log_amplitude = 0.0; end;
if ~isfield(parameter_dolphin,'X_log_amplitude'); parameter_dolphin.X_log_amplitude = 0.0; end;
if ~isfield(parameter_dolphin,'n_iteration_BtBn'); parameter_dolphin.n_iteration_BtBn = 32; end;
if ~isfield(parameter_dolphin,'MaxFunEvals_use_BtBn'); parameter_dolphin.MaxFunEvals_use_BtBn = 1024; end;
if ~isfield(parameter_dolphin,'MaxFunEvals_use_simultaneous'); parameter_dolphin.MaxFunEvals_use_simultaneous = 16; end;
if ~isfield(parameter_dolphin,'flag_regularize_eccentricity_simultaneous'); parameter_dolphin.flag_regularize_eccentricity_simultaneous = 1; end;
if ~isfield(parameter_dolphin,'rseed'); parameter_dolphin.rseed = 1; end;
%%%%%%%%;
dir_trunk = parameter_dolphin.dir_trunk;
dir_dolphin_mat = parameter_dolphin.dir_dolphin_mat;
tolerance_master = parameter_dolphin.tolerance_master;
flag_verbose = parameter_dolphin.flag_verbose;
flag_disp = parameter_dolphin.flag_disp;
dt_avg = parameter_dolphin.dt_avg;
n_i = parameter_dolphin.n_i;
n_j_factor = parameter_dolphin.n_j_factor;
ignore_factor = parameter_dolphin.ignore_factor;
T_ini = parameter_dolphin.T_ini;
T_max = parameter_dolphin.T_max;
B_log_amplitude = parameter_dolphin.B_log_amplitude;
C_log_amplitude = parameter_dolphin.C_log_amplitude;
X_log_amplitude = parameter_dolphin.X_log_amplitude;
n_iteration_BtBn = parameter_dolphin.n_iteration_BtBn;
MaxFunEvals_use_BtBn = parameter_dolphin.MaxFunEvals_use_BtBn;
MaxFunEvals_use_simultaneous = parameter_dolphin.MaxFunEvals_use_simultaneous;
flag_regularize_eccentricity_simultaneous = parameter_dolphin.flag_regularize_eccentricity_simultaneous;
rseed = parameter_dolphin.rseed;
%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [entering %s],',str_thisfunction)); end;

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

if ~exist(dir_dolphin_mat,'dir'); disp(sprintf(' %% mkdir %s',dir_dolphin_mat)); mkdir(dir_dolphin_mat); end;
fname_pre = sprintf('%s/test_dolphin_%s',dir_dolphin_mat,str_infix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% generate SDE data for each individual. ;
%%%%%%%%;
n_x = 2;
rng(rseed);
%%%%;
flag_stop = 0;
while ~flag_stop;
A_tru_xx__ = randn(n_x,n_x);
flag_stop = (trace(A_tru_xx__)<0) & (det(A_tru_xx__)>0); %<-- ensures eigenvalues of A_tru_xx__ have negative real part. ;
end;%while ~flag_stop;
A_tru_xx__ = A_tru_xx__/fnorm(A_tru_xx__);
%%%%;
n_a = 1; a_tru_xa__ = zeros(n_x,n_a);
B_tru_omega = 2*pi*rand();
B_tru_l0 = randn() + B_log_amplitude;
B_tru_l1 = randn() + B_log_amplitude;
C_tru_omega = 2*pi*rand();
C_tru_l0 = randn() + C_log_amplitude;
C_tru_l1 = randn() + C_log_amplitude;
X_ini_x_ = randn(n_x,1) .* exp(X_log_amplitude);
parameter_gen = struct('type','parameter');
parameter_gen.tolerance_master = tolerance_master;
parameter_gen.flag_verbose = flag_verbose;
parameter_gen.flag_disp = 0*flag_verbose;
parameter_gen.dt_avg = dt_avg;
parameter_gen.flag_discrete_vs_exponential = 0;
parameter_gen.T_ini = T_ini;
parameter_gen.T_max = T_max;
parameter_gen.n_j_factor = n_j_factor;
parameter_gen.ignore_factor = ignore_factor;
parameter_gen.rseed = rseed;
[ ...
,parameter_gen ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,X_tru_ixt___ ...
,BtBn_tru_xx__ ...
,BtBn_tru_inv_xx__ ...
,CtCn_tru_xx__ ...
,CtCn_tru_inv_xx__ ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,zlim_i2x___ ...
,P_tru_ixt___ ...
] = ...
SDE_generate_data_2_wrap_0( ...
 parameter_gen ...
,n_i ...
,n_x ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
,X_ini_x_ ...
);
%%%%%%%%;
% calculate average variation. ;
%%%%%%%%;
n_dt_sum = 0;
APdt_L2_sum = 0;
for ni=0:n_i-1;
n_t = n_t_i_(1+ni);
t_t_ = t_it__{1+ni};
dt_dt_ = diff(t_t_); n_dt = numel(dt_dt_);
P_tru_xt__ = P_tru_ixt___{1+ni};
APdt_xdt__ = bsxfun(@times,A_tru_xx__*P_tru_xt__(:,1:n_dt),reshape(dt_dt_,[1,n_dt]));
n_dt_sum = n_dt_sum + n_dt;
APdt_L2_sum = APdt_L2_sum + sum(APdt_xdt__.^2,'all');
clear n_t t_t_ dt_dt_ n_dt P_tru_xt__ APdt_xdt__;
end;%for ni=0:n_i-1;
APdt_L2_avg = APdt_L2_sum/n_dt_sum;
BidW_L2_avg = (exp(-B_tru_l0) + exp(-B_tru_l1))*dt_avg;
Cinv_L2_avg = (exp(-C_tru_l0) + exp(-C_tru_l1))*1;
%{
  %%%%;
  % Note average norm of BdW contribution is: ;
  %%%%;
  tmp0_xt__ = randn(n_x,1024*128); 
  tmp0 = mean(sum((BtBn_tru_inv_xx__*tmp0_xt__).*tmp0_xt__,1),2)*dt_avg;;
  tmp1_xt__ = sqrtm(BtBn_tru_inv_xx__)*tmp0_xt__*sqrt(dt_avg);
  tmp1 = mean(sum(tmp1_xt__.^2,1),2);
  tmp2 = (exp(-1*B_tru_l0)+exp(-1*B_tru_l1))*dt_avg;
  disp(sprintf(' %% tmp0 %0.6f tmp1 %0.6f tmp2 %0.6f',tmp0,tmp1,tmp2));
 %}
APdt_l2_avg = sqrt(APdt_L2_avg); BidW_l2_avg = sqrt(BidW_L2_avg); Cinv_l2_avg = sqrt(Cinv_L2_avg);
snr_A_vs_BC_avg = sqrt(APdt_L2_avg)./max(1e-12,sqrt(BidW_L2_avg + Cinv_L2_avg));
snr_AB_vs_C_avg = sqrt(APdt_L2_avg + BidW_L2_avg)./max(1e-12,sqrt(Cinv_L2_avg));
if (flag_verbose>0); disp(sprintf(' %% sum(n_j_i_) %.5d ; sum(n_t_i_) %.5d <-- snr_A_vs_BC_avg %0.6f snr_AB_vs_C_avg %0.6f (APdt_l2_avg %0.6f BidW_l2_avg %0.6f Cinv_l2_avg %0.6f)',sum(n_j_i_),sum(n_t_i_),snr_A_vs_BC_avg,snr_AB_vs_C_avg,APdt_l2_avg,BidW_l2_avg,Cinv_l2_avg)); end;

%%%%%%%%;
% recover all. ;
%%%%%%%%;
X_est_ixt___=[];
a_est_xa__=zeros(n_x,n_a);
A_est_xx__=zeros(n_x,n_x);
B_est_omega=0;
B_est_l0=-10;
B_est_l1=-10;
C_est_omega=0;
C_est_l0=-10;
C_est_l1=-10;
%%%%;
parameter_est = struct('type','parameter_est');
parameter_est.tolerance_master = tolerance_master;
parameter_est.flag_verbose = flag_verbose-1;
parameter_est.flag_disp = flag_verbose-1;
parameter_est.flag_regularize_eccentricity_BtBn = 1;
parameter_est.flag_regularize_eccentricity_simultaneous = flag_regularize_eccentricity_simultaneous;
parameter_est.n_iteration_BtBn = n_iteration_BtBn;
parameter_est.MaxFunEvals_use_BTBn = MaxFunEvals_use_BtBn;
parameter_est.MaxFunEvals_use_simultaneous = MaxFunEvals_use_simultaneous;
[ ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_est_ixt___ ...
,n_a ...
,a_est_xa__ ...
,A_est_xx__ ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
] = ...
SDE_nlp_ijXaABYC_update_all_1( ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_est_ixt___ ...
,n_a ...
,a_est_xa__ ...
,A_est_xx__ ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);
%%%%%%%%;
[ ...
 ~ ...
,BtBn_est_xx__ ...
] = ...
SDE_BtBn_0( ...
 [] ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
);
BtBn_est_inv_xx__ = pinv(BtBn_est_xx__,tolerance_master);
%%%%;
[ ...
 ~ ...
,CtCn_est_xx__ ...
] = ...
SDE_BtBn_0( ...
 [] ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);
CtCn_est_inv_xx__ = pinv(CtCn_est_xx__,tolerance_master);
%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
p_row = 2; p_col = n_i; np=0;
for prow=0:p_row-1;
for pcol=0:p_col-1;
ni = pcol;
t_t_ = t_it__{1+ni};
if prow==0; tmp_W_xt__ = X_tru_ixt___{1+ni}; tmp_str = 'X_tru_xt__'; end;
if prow==1; tmp_W_xt__ = X_est_ixt___{1+ni}; tmp_str = 'X_est_xt__'; end;
subplot(p_row,p_col,1+pcol+prow*p_col);
plot(t_t_,tmp_W_xt__(1,:),'r-',t_t_,tmp_W_xt__(2,:),'b-');
xlabel('time');ylabel(tmp_str,'Interpreter','none'); xlim([T_ini,T_max]);
%s = surfline_0(tmp_W_xt__(1,:),tmp_W_xt__(2,:),t_t_); set(s,'LineWidth',3);
%xlabel(sprintf('%s(1+0,:)',tmp_str),'Interpreter','none');
%ylabel(sprintf('%s(1+1,:)',tmp_str),'Interpreter','none');
title(sprintf('ni %d',ni),'Interpreter','none');
end;%for pcol=0:p_col-1;
end;%for prow=0:p_row-1;
sgtitle(sprintf('all: X_est_ixt___'),'Interpreter','none');
end;%if flag_disp>0;
%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
p_row = 2; p_col = n_i; np=0;
for prow=0:p_row-1;
for pcol=0:p_col-1;
ni = pcol;
t_t_ = t_it__{1+ni};
if prow==0; tmp_W_xt__ = X_tru_ixt___{1+ni}; tmp_str = 'X_tru_xt__'; end;
if prow==1; tmp_W_xt__ = X_est_ixt___{1+ni}; tmp_str = 'X_est_xt__'; end;
subplot(p_row,p_col,1+pcol+prow*p_col);
%plot(t_t_,tmp_W_xt__(1,:),'r-',t_t_,tmp_W_xt__(2,:),'b-');
%xlabel('time');ylabel(tmp_str,'Interpreter','none'); xlim([T_ini,T_max]);
s = surfline_0(tmp_W_xt__(1,:),tmp_W_xt__(2,:),t_t_); set(s,'LineWidth',3);
xlabel(sprintf('%s(1+0,:)',tmp_str),'Interpreter','none');
ylabel(sprintf('%s(1+1,:)',tmp_str),'Interpreter','none');
title(sprintf('ni %d',ni),'Interpreter','none');
end;%for pcol=0:p_col-1;
end;%for prow=0:p_row-1;
sgtitle(sprintf('all: X_est_ixt___'),'Interpreter','none');
end;%if flag_disp>0;
%%%%%%%%;

rfnorm_A = fnorm(A_tru_xx__ - A_est_xx__)/max(1e-12,fnorm(A_tru_xx__));
corr_A = corr(A_tru_xx__(:),A_est_xx__(:));
rfnorm_B = fnorm(BtBn_tru_xx__ - BtBn_est_xx__)/max(1e-12,fnorm(BtBn_tru_xx__));
corr_B = corr(BtBn_tru_xx__(:),BtBn_est_xx__(:));
rfnorm_B_inv = fnorm(BtBn_tru_inv_xx__ - BtBn_est_inv_xx__)/max(1e-12,fnorm(BtBn_tru_inv_xx__));
corr_B_inv = corr(BtBn_tru_inv_xx__(:),BtBn_est_inv_xx__(:));
rfnorm_C = fnorm(CtCn_tru_xx__ - CtCn_est_xx__)/max(1e-12,fnorm(CtCn_tru_xx__));
corr_C = corr(CtCn_tru_xx__(:),CtCn_est_xx__(:));
rfnorm_C_inv = fnorm(CtCn_tru_inv_xx__ - CtCn_est_inv_xx__)/max(1e-12,fnorm(CtCn_tru_inv_xx__));
corr_C_inv = corr(CtCn_tru_inv_xx__(:),CtCn_est_inv_xx__(:));

if (flag_verbose>0);
disp(sprintf(' %% A_tru_xx__ vs A_est_xx__: %+0.6f',rfnorm_A));
disp(sprintf(' %% corr(A_tru_xx__(:),A_est_xx__(:)): %+0.6f',corr_A));
disp(sprintf(' %% BtBn_tru_xx__ vs BtBn_est_xx__: %+0.6f',rfnorm_B));
disp(sprintf(' %% corr(BtBn_tru_xx__(:),BtBn_est_xx__(:)): %+0.6f',corr_B));
disp(sprintf(' %% BtBn_tru_inv_xx__ vs BtBn_est_inv_xx__: %+0.6f',rfnorm_B_inv));
disp(sprintf(' %% corr(BtBn_tru_inv_xx__(:),BtBn_est_inv_xx__(:)): %+0.6f',corr_B_inv));
disp(sprintf(' %% CtCn_tru_xx__ vs CtCn_est_xx__: %+0.6f',rfnorm_C));
disp(sprintf(' %% corr(CtCn_tru_xx__(:),CtCn_est_xx__(:)): %+0.6f',corr_C));
disp(sprintf(' %% CtCn_tru_inv_xx__ vs CtCn_est_inv_xx__: %+0.6f',rfnorm_C_inv));
disp(sprintf(' %% corr(CtCn_tru_inv_xx__(:),CtCn_est_inv_xx__(:)): %+0.6f',corr_C_inv));
end;%if (flag_verbose>0); 

save(fname_mat ...
,'parameter_dolphin' ...
,'dt_avg' ...
,'n_j_factor' ...
,'ignore_factor' ...
,'T_ini' ...
,'T_max' ...
,'B_log_amplitude' ...
,'C_log_amplitude' ...
,'X_log_amplitude' ...
,'MaxFunEvals_use_simultaneous' ...
,'flag_regularize_eccentricity_simultaneous' ...
,'rseed' ...
,'str_infix' ...
,'n_x' ...
,'A_tru_xx__' ...
,'n_a' ...
,'a_tru_xa__' ...
,'B_tru_omega' ...
,'B_tru_l0' ...
,'B_tru_l1' ...
,'C_tru_omega' ...
,'C_tru_l0' ...
,'C_tru_l1' ...
,'X_ini_x_' ...
,'parameter_gen' ...
,'n_t_i_' ...
,'BtBn_tru_xx__' ...
,'BtBn_tru_inv_xx__' ...
,'CtCn_tru_xx__' ...
,'CtCn_tru_inv_xx__' ...
,'n_j_i_' ...
,'APdt_L2_avg' ...
,'BidW_L2_avg' ...
,'Cinv_L2_avg' ...
,'APdt_l2_avg' ...
,'BidW_l2_avg' ...
,'Cinv_l2_avg' ...
,'snr_A_vs_BC_avg' ...
,'snr_AB_vs_C_avg' ...
,'a_est_xa__' ...
,'A_est_xx__' ...
,'B_est_omega' ...
,'B_est_l0' ...
,'B_est_l1' ...
,'C_est_omega' ...
,'C_est_l0' ...
,'C_est_l1' ...
,'parameter_est' ...
,'BtBn_est_xx__' ...
,'BtBn_est_inv_xx__' ...
,'CtCn_est_xx__' ...
,'CtCn_est_inv_xx__' ...
,'rfnorm_A' ...
,'corr_A' ...
,'rfnorm_B' ...
,'corr_B' ...
,'rfnorm_B_inv' ...
,'corr_B_inv' ...
,'rfnorm_C' ...
,'corr_C' ...
,'rfnorm_C_inv' ...
,'corr_C_inv' ...
);
close_fname_tmp(fname_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s],',str_thisfunction)); end;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;


