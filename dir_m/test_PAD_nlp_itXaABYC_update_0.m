%%%%%%%%;
tolerance_master = 1e-6;
flag_verbose = 1;
flag_disp = flag_verbose; nf=0;
rseed = 1;
%%%%%%%%;

%%%%%%%%;
% generate SDE data for each individual. ;
%%%%%%%%;
n_x = 2; T_ini = 0; T_max = 256; dt_0in = 1/16;
A_tru_xx__ = [-0.1 , +1.0   ; ...
              -1.0 , -0.2 ] ;
n_a = 3;
a_tru_xa__ = [ +1.00 , -2.00/T_max , +9.50/T_max.^2   ; ...
               -0.50 , +1.25/T_max , -9.25/T_max.^2 ] ;
B_tru_omega = +pi/3;
B_tru_l0 = +1.0;
B_tru_l1 = -0.5;
[ ...
 ~ ...
,BtBn_tru_xx__ ...
] = ...
PAD_BtBn_0( ...
 [] ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
);
C_tru_omega = -2*pi/5;
C_tru_l0 = +0.3;
C_tru_l1 = -0.8;
[ ...
 ~ ...
,CtCn_tru_xx__ ...
] = ...
PAD_BtBn_0( ...
 [] ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
);
BtBn_tru_inv_xx__ = pinv(BtBn_tru_xx__,tolerance_master);
CtCn_tru_inv_xx__ = pinv(CtCn_tru_xx__,tolerance_master);
%%%%%%%%;
n_i = 3; %<-- number of individuals. ;
n_t_i_ = zeros(n_i,1);
t_it__ = cell(n_i,1);
Q_tru_ixt___ = cell(n_i,1);
P_tru_ixt___ = cell(n_i,1);
X_tru_ixt___ = cell(n_i,1);
Y_tru_ixt___ = cell(n_i,1);
ignore_Y_tru_ixt___ = cell(n_i,1);
%%%%;
parameter_gen = struct('type','parameter');
parameter_gen.dt_avg = dt_0in;
parameter_gen.flag_discrete_vs_exponential = 1;
parameter_gen.T_ini = T_ini;
parameter_gen.T_max = T_max;
for ni=0:n_i-1;
parameter_gen.rseed = 1+ni;
[ ...
 parameter_gen ...
,n_t ...
,t_t_ ...
,Q_tru_xt__ ...
,P_tru_xt__ ...
,X_tru_xt__ ...
,Y_tru_xt__ ...
,R_avg ...
] = ...
SDE_generate_data_1( ...
 parameter_gen ...
,n_x ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,BtBn_tru_xx__ ...
,BtBn_tru_inv_xx__ ...
,CtCn_tru_xx__ ...
,CtCn_tru_inv_xx__ ...
);
%%%%;
if flag_disp>1;
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
p_row = 2; p_col = 4; np=0;
for pcol=0:p_col-1;
if pcol==0; tmp_W_xt__ = Q_tru_xt__; tmp_str = 'Q_tru_xt__'; end;
if pcol==1; tmp_W_xt__ = P_tru_xt__; tmp_str = 'P_tru_xt__'; end;
if pcol==2; tmp_W_xt__ = X_tru_xt__; tmp_str = 'X_tru_xt__'; end;
if pcol==3; tmp_W_xt__ = Y_tru_xt__; tmp_str = 'Y_tru_xt__'; end;
subplot(p_row,p_col,1+pcol+0*p_col);
plot(t_t_,tmp_W_xt__(1,:),'r-',t_t_,tmp_W_xt__(2,:),'b-');
xlabel('time');ylabel(tmp_str,'Interpreter','none'); xlim([T_ini,T_max]);
subplot(p_row,p_col,1+pcol+1*p_col);
s = surfline_0(tmp_W_xt__(1,:),tmp_W_xt__(2,:),t_t_); set(s,'LineWidth',3);
xlabel(sprintf('%s(1+0,:)',tmp_str),'Interpreter','none');
ylabel(sprintf('%s(1+1,:)',tmp_str),'Interpreter','none');
end;%for pcol=0:p_col-1;
sgtitle(sprintf('ni %d',ni),'Interpreter','none');
end;%if flag_disp;
%%%%;
n_t_i_(1+ni) = n_t;
t_it__{1+ni} = t_t_;
Q_tru_ixt___{1+ni} = Q_tru_xt__;
P_tru_ixt___{1+ni} = P_tru_xt__;
X_tru_ixt___{1+ni} = X_tru_xt__;
Y_tru_ixt___{1+ni} = Y_tru_xt__;
rng(parameter_gen.rseed); ignore_Y_tru_ixt___{1+ni} = rand(n_x,n_t)<0.5;
clear n_t t_t_ Q_tru_xt__ P_tru_xt__ X_tru_xt__ Y_tru_xt__ R_avg ;
%%%%;
end;%for ni=0:n_i-1;

%%%%%%%%;
% Recover CtCn_xx__. ;
%%%%%%%%;
parameter_est = struct('type','parameter_est');
parameter_est.tolerance_master = tolerance_master;
parameter_est.flag_verbose = flag_verbose;
parameter_est.str_update = 'CtCn_xx__';
[ ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_tru_ixt____ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
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
,X_tru_ixt___ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,ignore_Y_tru_ixt___ ...
,Y_tru_ixt___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
);
%%%%;
[ ...
 ~ ...
,CtCn_est_xx__ ...
] = ...
PAD_BtBn_0( ...
 [] ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);
%%%%;
if flag_verbose;
disp(sprintf(' %% nlp_tXYC_pre %0.6f --> nlp_tXYC_pos %0.6f',parameter_est.nlp_tXYC_pre,parameter_est.nlp_tXYC_pos));
disp(sprintf(' %% CtCn_tru_xx__: ')); disp(CtCn_tru_xx__);
disp(sprintf(' %% CtCn_est_xx__: ')); disp(CtCn_est_xx__);
end;%if flag_verbose;
%%%%%%%%;

%%%%%%%%;
% Recover BtBn_xx__. ;
%%%%%%%%;
parameter_est = struct('type','parameter_est');
parameter_est.tolerance_master = tolerance_master;
parameter_est.flag_verbose = flag_verbose;
parameter_est.str_update = 'BtBn_xx__';
[ ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_tru_ixt____ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
,ignore_Y_tru_ixt____ ...
,Y_tru_ixt____ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
] = ...
PAD_nlp_itXaABYC_update_0( ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_tru_ixt___ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,ignore_Y_tru_ixt___ ...
,Y_tru_ixt___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
);
%%%%;
[ ...
 ~ ...
,BtBn_est_xx__ ...
] = ...
PAD_BtBn_0( ...
 [] ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
);
%%%%;
if flag_verbose;
disp(sprintf(' %% nlp_dtZPB_pre %0.6f --> nlp_dtZPB_pos %0.6f',parameter_est.nlp_dtZPB_pre,parameter_est.nlp_dtZPB_pos));
disp(sprintf(' %% BtBn_tru_xx__: ')); disp(BtBn_tru_xx__);
disp(sprintf(' %% BtBn_est_xx__: ')); disp(BtBn_est_xx__);
end;%if flag_verbose;
%%%%%%%%;

%%%%%%%%;
% Recover a_xa__ and A_xx__. ;
%%%%%%%%;
parameter_est = struct('type','parameter_est');
parameter_est.tolerance_master = tolerance_master;
parameter_est.flag_verbose = flag_verbose-1;
parameter_est.str_update = 'A_xx__';
[ ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_tru_ixt____ ...
,n_a ...
,a_est_xa__ ...
,A_est_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,ignore_Y_tru_ixt____ ...
,Y_tru_ixt____ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
] = ...
PAD_nlp_itXaABYC_update_0( ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_tru_ixt___ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,ignore_Y_tru_ixt___ ...
,Y_tru_ixt___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
);
%%%%;
if flag_verbose;
disp(sprintf(' %% nlp_dtXaAB_sum_pre %0.6f --> nlp_dtXaAB_sum_pos %0.6f',parameter_est.nlp_dtXaAB_sum_pre,parameter_est.nlp_dtXaAB_sum_pos));
disp(sprintf(' %% a_tru_xa__: ')); disp(a_tru_xa__);
disp(sprintf(' %% a_est_xa__: ')); disp(a_est_xa__);
disp(sprintf(' %% A_tru_xx__: ')); disp(A_tru_xx__);
disp(sprintf(' %% A_est_xx__: ')); disp(A_est_xx__);
end;%if flag_verbose;
%%%%%%%%;

%%%%%%%%;
% Recover X_xt__. ;
%%%%%%%%;
parameter_est = struct('type','parameter_est');
parameter_est.tolerance_master = tolerance_master;
parameter_est.flag_verbose = flag_verbose-1;
parameter_est.str_update = 'X_xt__';
[ ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_est_ixt___ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,ignore_Y_tru_ixt____ ...
,Y_tru_ixt____ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
] = ...
PAD_nlp_itXaABYC_update_0( ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_tru_ixt___ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,ignore_Y_tru_ixt___ ...
,Y_tru_ixt___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
);
%%%%;
if flag_verbose;
disp(sprintf(' %% nlp_tXaABYC_sum_pre %0.6f --> nlp_tXaABYC_sum_pos %0.6f',parameter_est.nlp_tXaABYC_sum_pre,parameter_est.nlp_tXaABYC_sum_pos));
disp(sprintf(' %% nlp_tXaABYC_integrated_sum_pre %0.6f --> nlp_tXaABYC_integrated_sum_pos %0.6f',parameter_est.nlp_tXaABYC_integrated_sum_pre,parameter_est.nlp_tXaABYC_integrated_sum_pos));
end;%if flag_verbose;
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
sgtitle(sprintf('X_est_ixt___'),'Interpreter','none');
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
sgtitle(sprintf('X_est_ixt___'),'Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% recover all. ;
%%%%%%%%;
X_est_ixt___ = X_tru_ixt___;
a_est_xa__ = a_tru_xa__;
A_est_xx__ = A_tru_xx__;
B_est_omega = B_tru_omega;
B_est_l0 = B_tru_l0;
B_est_l1 = B_tru_l1;
C_est_omega = C_tru_omega;
C_est_l0 = C_tru_l0;
C_est_l1 = C_tru_l1;
%%%%;
X_est_ixt___=[];
a_est_xa__=[];
A_est_xx__=[];
B_est_omega=[];
B_est_l0=[];
B_est_l1=[];
C_est_omega=[];
C_est_l0=[];
C_est_l1=[];
%%%%;
parameter_est = struct('type','parameter_est');
parameter_est.tolerance_master = tolerance_master;
parameter_est.flag_verbose = flag_verbose-1;
parameter_est.str_update = 'CtCn_xx__ BtBn_xx__ a_xa__ A_xx__ X_xt__';
parameter_est.flag_regularize_eccentricity = 1;
parameter_est.MaxFunEvals_use = 1024;
n_iteration = 6;
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
parameter_est.flag_regularize_eccentricity = 0;
parameter_est.MaxFunEvals_use = 1024*4;
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
disp(sprintf(' %% nlp_itaABYC_pre %0.6f --> nlp_itaABYC_pos %0.6f',parameter_est.nlp_itaABYC_pre,parameter_est.nlp_itaABYC_pos));
disp(sprintf(' %% nlp_itaABYC_pre %0.6f --> nlp_itaABYC_pos %0.6f',parameter_est.nlp_itaABYC_pre,parameter_est.nlp_itaABYC_pos));
end;%if flag_verbose;
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


