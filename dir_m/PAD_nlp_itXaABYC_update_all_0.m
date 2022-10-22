function ...
[ ...
 parameter ...
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
,ignore_Y_tru_ixt___ ...
,Y_tru_ixt___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
] = ...
PAD_nlp_itXaABYC_update_all_0( ...
 parameter ...
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
,ignore_Y_tru_ixt___ ...
,Y_tru_ixt___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);

str_thisfunction = 'PAD_nlp_itXaABYC_update_all_0';
nf=0;

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
disp(sprintf(' %% see: test_%s_0.m',str_thisfunction));
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_i=[]; end; na=na+1;
if (nargin<1+na); n_t_i_=[]; end; na=na+1;
if (nargin<1+na); t_it__=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); X_est_ixt___=[]; end; na=na+1;
if (nargin<1+na); n_a=[]; end; na=na+1;
if (nargin<1+na); a_est_xa__=[]; end; na=na+1;
if (nargin<1+na); A_est_xx__=[]; end; na=na+1;
if (nargin<1+na); B_est_omega=[]; end; na=na+1;
if (nargin<1+na); B_est_l0=[]; end; na=na+1;
if (nargin<1+na); B_est_l1=[]; end; na=na+1;
if (nargin<1+na); ignore_Y_tru_ixt___=[]; end; na=na+1;
if (nargin<1+na); Y_tru_ixt___=[]; end; na=na+1;
if (nargin<1+na); C_est_omega=[]; end; na=na+1;
if (nargin<1+na); C_est_l0=[]; end; na=na+1;
if (nargin<1+na); C_est_l1=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
if ~isfield(parameter,'flag_regularize_eccentricity_BtBn'); parameter.flag_regularize_eccentricity_BtBn = 1; end;
if ~isfield(parameter,'flag_regularize_eccentricity_simultaneous'); parameter.flag_regularize_eccentricity_simultaneous = 0; end;
if ~isfield(parameter,'n_iteration_BtBn'); parameter.n_iteration_BtBn = 16; end;
if ~isfield(parameter,'MaxFunEvals_use_BtBn'); parameter.MaxFunEvals_use_BtBn = 1024; end;
if ~isfield(parameter,'MaxFunEvals_use_simultaneous'); parameter.MaxFunEvals_use_simultaneous = 1024; end;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;
flag_disp = parameter.flag_disp;
flag_regularize_eccentricity_BtBn = parameter.flag_regularize_eccentricity_BtBn;
flag_regularize_eccentricity_simultaneous = parameter.flag_regularize_eccentricity_simultaneous;
n_iteration_BtBn = parameter.n_iteration_BtBn;
MaxFunEvals_use_BtBn = parameter.MaxFunEvals_use_BtBn;
MaxFunEvals_use_simultaneous = parameter.MaxFunEvals_use_simultaneous;

%%%%%%%%;
% recover all. ;
%%%%%%%%;
parameter_est = struct('type','parameter_est');
parameter_est.tolerance_master = tolerance_master;
parameter_est.flag_verbose = flag_verbose-1;
parameter_est.str_update = 'CtCn_xx__ BtBn_xx__ a_xa__ A_xx__ X_xt__';
parameter_est.flag_regularize_eccentricity = flag_regularize_eccentricity_BtBn;
parameter_est.MaxFunEvals_use = MaxFunEvals_use_BtBn;
n_iteration = n_iteration_BtBn;
nlp_tXYC_l_ = zeros(n_iteration,1);
nlp_dtZPB_l_ = zeros(n_iteration,1);
nlp_dtXaAB_sum_l_ = zeros(n_iteration,1);
nlp_tXaABYC_sum_l_ = zeros(n_iteration,1);
nlp_tXaABYC_integrated_sum_l_ = zeros(n_iteration,1);
%%%%;
flag_continue = 1; niteration=0; nlp_tXaABYC_integrated_sum_pos=+Inf;
while flag_continue;
nlp_tXaABYC_integrated_sum_pre = nlp_tXaABYC_integrated_sum_pos;
[ ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_tmp_ixt___ ...
,n_a ...
,a_tmp_xa__ ...
,A_tmp_xx__ ...
,B_tmp_omega ...
,B_tmp_l0 ...
,B_tmp_l1 ...
,ignore_Y_tru_ixt____ ...
,Y_tru_ixt____ ...
,C_tmp_omega ...
,C_tmp_l0 ...
,C_tmp_l1 ...
] = ...
PAD_nlp_itXaABYC_update_0( ...
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
,ignore_Y_tru_ixt___ ...
,Y_tru_ixt___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);
%%%%;
if flag_verbose;
disp(sprintf(' %% processing niteration %d/%d',niteration,n_iteration));
[~,BtBn_est_xx__] =PAD_BtBn_0([],B_est_omega,B_est_l0,B_est_l1);
[~,CtCn_est_xx__] =PAD_BtBn_0([],C_est_omega,C_est_l0,C_est_l1);
disp(sprintf(' %% B_omega %+0.2f, B_l0 %+0.2f, B_l1 %+0.2f',B_est_omega,B_est_l0,B_est_l1));
disp(sprintf(' %% C_omega %+0.2f, C_l0 %+0.2f, C_l1 %+0.2f',C_est_omega,C_est_l0,C_est_l1));
if isfield(parameter_est,'nlp_tXYC_pos'); disp(sprintf(' %% niteration %d/%d: nlp_tXYC_pre %0.6f --> nlp_tXYC_pos %0.6f',niteration,n_iteration,parameter_est.nlp_tXYC_pre,parameter_est.nlp_tXYC_pos)); nlp_tXYC_l_(1+niteration) = parameter_est.nlp_tXYC_pos; end;
if isfield(parameter_est,'nlp_dtZPB_pos'); disp(sprintf(' %% niteration %d/%d: nlp_dtZPB_pre %0.6f --> nlp_dtZPB_pos %0.6f',niteration,n_iteration,parameter_est.nlp_dtZPB_pre,parameter_est.nlp_dtZPB_pos)); nlp_dtZPB_l_(1+niteration) = parameter_est.nlp_dtZPB_pos; end;
if isfield(parameter_est,'nlp_dtXaAB_sum_pos'); disp(sprintf(' %% niteration %d/%d: nlp_dtXaAB_sum_pre %0.6f --> nlp_dtXaAB_sum_pos %0.6f',niteration,n_iteration,parameter_est.nlp_dtXaAB_sum_pre,parameter_est.nlp_dtXaAB_sum_pos)); nlp_dtXaAB_sum_l_(1+niteration) = parameter_est.nlp_dtXaAB_sum_pos; end;
if isfield(parameter_est,'nlp_tXaABYC_sum_pos'); disp(sprintf(' %% niteration %d/%d: nlp_tXaABYC_sum_pre %0.6f --> nlp_tXaABYC_sum_pos %0.6f',niteration,n_iteration,parameter_est.nlp_tXaABYC_sum_pre,parameter_est.nlp_tXaABYC_sum_pos)); nlp_tXaABYC_sum_l_(1+niteration) = parameter_est.nlp_tXaABYC_sum_pos; end;
if isfield(parameter_est,'nlp_tXaABYC_integrated_sum_pos'); disp(sprintf(' %% niteration %d/%d: nlp_tXaABYC_integrated_sum_pre %0.6f --> nlp_tXaABYC_integrated_sum_pos %0.6f',niteration,n_iteration,parameter_est.nlp_tXaABYC_integrated_sum_pre,parameter_est.nlp_tXaABYC_integrated_sum_pos)); nlp_tXaABYC_integrated_sum_l_(1+niteration) = parameter_est.nlp_tXaABYC_integrated_sum_pos; end;
end;%if flag_verbose;
%%%%;
flag_continue = 0;
nlp_tXaABYC_integrated_sum_pos = parameter_est.nlp_tXaABYC_integrated_sum_pos;
if (nlp_tXaABYC_integrated_sum_pos >= nlp_tXaABYC_integrated_sum_pre);
if (flag_verbose);
disp(sprintf(' %% niteration %d/%d: nlp_tXaABYC_integrated_sum_pos %0.6f >= nlp_tXaABYC_integrated_sum_pre %0.6f',niteration,n_iteration,nlp_tXaABYC_integrated_sum_pos,nlp_tXaABYC_integrated_sum_pre));
end;%if (flag_verbose);
flag_continue = 0; %<-- stopping. ;
end;%if (nlp_tXaABYC_integrated_sum_pos >= nlp_tXaABYC_integrated_sum_pre);
if (nlp_tXaABYC_integrated_sum_pos <  nlp_tXaABYC_integrated_sum_pre);
if (flag_verbose);
disp(sprintf(' %% niteration %d/%d: nlp_tXaABYC_integrated_sum_pos %0.6f <  nlp_tXaABYC_integrated_sum_pre %0.6f',niteration,n_iteration,nlp_tXaABYC_integrated_sum_pos,nlp_tXaABYC_integrated_sum_pre));
end;%if (flag_verbose);
X_est_ixt___ = X_tmp_ixt___;
a_est_xa__ = a_tmp_xa__;
A_est_xx__ = A_tmp_xx__;
B_est_omega = B_tmp_omega;
B_est_l0 = B_tmp_l0;
B_est_l1 = B_tmp_l1;
C_est_omega = C_tmp_omega;
C_est_l0 = C_tmp_l0;
C_est_l1 = C_tmp_l1;
flag_continue = 1; %<-- continue. ;
end;%if (nlp_tXaABYC_integrated_sum_pos <  nlp_tXaABYC_integrated_sum_pre);
flag_continue = flag_continue & (niteration<n_iteration) ;
niteration = niteration + 1;
end;%while flag_continue;
%%%%;
if flag_verbose;
disp(sprintf(' %% stopping at niteration %d/%d',niteration,n_iteration));
[~,BtBn_est_xx__] =PAD_BtBn_0([],B_est_omega,B_est_l0,B_est_l1);
[~,CtCn_est_xx__] =PAD_BtBn_0([],C_est_omega,C_est_l0,C_est_l1);
disp(sprintf(' %% B_omega %+0.2f, B_l0 %+0.2f, B_l1 %+0.2f',B_est_omega,B_est_l0,B_est_l1));
disp(sprintf(' %% C_omega %+0.2f, C_l0 %+0.2f, C_l1 %+0.2f',C_est_omega,C_est_l0,C_est_l1));
if isfield(parameter_est,'nlp_tXYC_pos'); disp(sprintf(' %% niteration %d/%d: nlp_tXYC_pre %0.6f --> nlp_tXYC_pos %0.6f',niteration,n_iteration,parameter_est.nlp_tXYC_pre,parameter_est.nlp_tXYC_pos)); nlp_tXYC_l_(1+niteration) = parameter_est.nlp_tXYC_pos; end;
if isfield(parameter_est,'nlp_dtZPB_pos'); disp(sprintf(' %% niteration %d/%d: nlp_dtZPB_pre %0.6f --> nlp_dtZPB_pos %0.6f',niteration,n_iteration,parameter_est.nlp_dtZPB_pre,parameter_est.nlp_dtZPB_pos)); nlp_dtZPB_l_(1+niteration) = parameter_est.nlp_dtZPB_pos; end;
if isfield(parameter_est,'nlp_dtXaAB_sum_pos'); disp(sprintf(' %% niteration %d/%d: nlp_dtXaAB_sum_pre %0.6f --> nlp_dtXaAB_sum_pos %0.6f',niteration,n_iteration,parameter_est.nlp_dtXaAB_sum_pre,parameter_est.nlp_dtXaAB_sum_pos)); nlp_dtXaAB_sum_l_(1+niteration) = parameter_est.nlp_dtXaAB_sum_pos; end;
if isfield(parameter_est,'nlp_tXaABYC_sum_pos'); disp(sprintf(' %% niteration %d/%d: nlp_tXaABYC_sum_pre %0.6f --> nlp_tXaABYC_sum_pos %0.6f',niteration,n_iteration,parameter_est.nlp_tXaABYC_sum_pre,parameter_est.nlp_tXaABYC_sum_pos)); nlp_tXaABYC_sum_l_(1+niteration) = parameter_est.nlp_tXaABYC_sum_pos; end;
if isfield(parameter_est,'nlp_tXaABYC_integrated_sum_pos'); disp(sprintf(' %% niteration %d/%d: nlp_tXaABYC_integrated_sum_pre %0.6f --> nlp_tXaABYC_integrated_sum_pos %0.6f',niteration,n_iteration,parameter_est.nlp_tXaABYC_integrated_sum_pre,parameter_est.nlp_tXaABYC_integrated_sum_pos)); nlp_tXaABYC_integrated_sum_l_(1+niteration) = parameter_est.nlp_tXaABYC_integrated_sum_pos; end;
end;%if flag_verbose;
%%%%;
if flag_verbose; disp(sprintf(' %% polishing with simultaneous iteration: ')); end;
parameter_est.tolerance_master = tolerance_master;
parameter_est.flag_verbose = flag_verbose-0;
parameter_est.str_update = 'simultaneous';
parameter_est.flag_regularize_eccentricity = flag_regularize_eccentricity_simultaneous;
parameter_est.MaxFunEvals_use = MaxFunEvals_use_simultaneous;
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
,ignore_Y_tru_ixt____ ...
,Y_tru_ixt____ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
] = ...
PAD_nlp_itXaABYC_update_0( ...
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
,ignore_Y_tru_ixt___ ...
,Y_tru_ixt___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);
%%%%;
if flag_verbose;
if isfield(parameter_est,'nlp_itaABYC_pre') & isfield(parameter_est,'nlp_itaABYC_pos');
disp(sprintf(' %% nlp_itaABYC_pre %0.6f --> nlp_itaABYC_pos %0.6f',parameter_est.nlp_itaABYC_pre,parameter_est.nlp_itaABYC_pos));
end;%if isfield(parameter_est,'nlp_itaABYC_pre') & isfield(parameter_est,'nlp_itaABYC_pos');
end;%if flag_verbose;
%%%%;
if isfield(parameter_est,'nlp_itaABYC_pos');
parameter.nlp_itaABYC = parameter_est.nlp_itaABYC_pos;
end;%if isfield(parameter_est,'nlp_itaABYC_pos');
%%%%;

%%%%;
if flag_disp;
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
end;%if flag_disp;
%%%%;
if flag_disp;
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
end;%if flag_disp;
%%%%%%%%;
