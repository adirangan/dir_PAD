%%%%%%%%;
% Built from test_SDE_nlp_ijXaABYC_update_1.m. ;
% Using statistics from dolphin data. ;
%%%%%%%%;

%%%%%%%%;
flag_verbose = 1;
flag_disp = 0*flag_verbose; nf=0;
tolerance_master = 1e-6;
dolphin_tolerance_master = tolerance_master;
dolphin_dt_avg = 0.257; %<-- matches d00 data. ;
dolphin_n_i = 1;%dolphin_n_i = 144; %<-- matches d00 data. ;
dolphin_n_j_factor = 1.02;
dolphin_ignore_factor = 0.02;
dolphin_T_ini =  0;
dolphin_T_max = 55*27*1;%dolphin_T_max = 55*144;%dolphin_T_max = 55;
dolphin_MaxFunEvals_use_simultaneous = 0;
dolphin_rseed = 2;
%%%%%%%%;

%%%%%%%%;
% generate SDE data for each individual. ;
%%%%%%%%;
n_x = 2;
rng(dolphin_rseed);
A_tru_xx__ = randn(n_x,n_x); [VA__,DA__] = eig(A_tru_xx__);
DA_ = diag(DA__); DA_ = min(-0.1,real(DA_)) + i*imag(DA_); DA__ = diag(DA_);
A_tru_xx__ = real(VA__*DA__*pinv(VA__,dolphin_tolerance_master));
n_a = 1; a_tru_xa__ = zeros(n_x,n_a);
B_tru_omega = 2*pi*rand();
B_tru_l0 = randn();
B_tru_l1 = randn();
C_tru_omega = 2*pi*rand();
C_tru_l0 = randn();
C_tru_l1 = randn();
parameter_gen = struct('type','parameter');
parameter_gen.tolerance_master = dolphin_tolerance_master;
parameter_gen.flag_verbose = flag_verbose;
parameter_gen.dt_avg = dolphin_dt_avg;
parameter_gen.flag_discrete_vs_exponential = 0;
parameter_gen.T_ini = dolphin_T_ini;
parameter_gen.T_max = dolphin_T_max;
parameter_gen.n_j_factor = dolphin_n_j_factor;
parameter_gen.ignore_factor = dolphin_ignore_factor;
parameter_gen.rseed = dolphin_rseed;
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
] = ...
SDE_generate_data_2_wrap_0( ...
 parameter_gen ...
,dolphin_n_i ...
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
);


%%%%%%%%;
% Recover CtCn_xx__. ;
%%%%%%%%;
C_est_omega = 0*C_tru_omega;
C_est_l0 = 0*C_tru_l0;
C_est_l1 = 0*C_tru_l1;
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
,X_tru_ixt___ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
] = ...
SDE_nlp_ijXaABYC_update_1( ...
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
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);
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
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% nlp_ijXYC_pre %0.6f --> nlp_ijXYC_pos %0.6f',parameter_est.nlp_ijXYC_pre,parameter_est.nlp_ijXYC_pos));
disp(sprintf(' %% CtCn_tru_xx__: ')); disp(CtCn_tru_xx__);
disp(sprintf(' %% CtCn_est_xx__: ')); disp(CtCn_est_xx__);
end;%if (flag_verbose>0);
%%%%%%%%;

%%%%%%%%;
% Recover BtBn_xx__. ;
%%%%%%%%;
B_est_omega = 0*B_tru_omega;
B_est_l0 = 0*B_tru_l0;
B_est_l1 = 0*B_tru_l1;
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
,X_tru_ixt___ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
] = ...
SDE_nlp_ijXaABYC_update_1( ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_tru_ixt___ ...
,n_a ...
,a_tru_xa__ ...
,A_tru_xx__ ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
);
%%%%;
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
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% nlp_idtZPB_pre %0.6f --> nlp_idtZPB_pos %0.6f',parameter_est.nlp_idtZPB_pre,parameter_est.nlp_idtZPB_pos));
disp(sprintf(' %% BtBn_tru_xx__: ')); disp(BtBn_tru_xx__);
disp(sprintf(' %% BtBn_est_xx__: ')); disp(BtBn_est_xx__);
end;%if (flag_verbose>0);
%%%%%%%%;

%%%%%%%%;
% Recover a_xa__ and A_xx__. ;
%%%%%%%%;
a_est_xa__ = zeros(n_x,n_a);
A_est_xx__ = zeros(n_x,n_x);
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
,X_tru_ixt___ ...
,n_a ...
,a_est_xa__ ...
,A_est_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
] = ...
SDE_nlp_ijXaABYC_update_1( ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_tru_ixt___ ...
,n_a ...
,a_est_xa__ ...
,A_est_xx__ ...
,B_tru_omega ...
,B_tru_l0 ...
,B_tru_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
);
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% nlp_idtXaAB_pre %0.6f --> nlp_idtXaAB_pos %0.6f',parameter_est.nlp_idtXaAB_pre,parameter_est.nlp_idtXaAB_pos));
disp(sprintf(' %% a_tru_xa__: ')); disp(a_tru_xa__);
disp(sprintf(' %% a_est_xa__: ')); disp(a_est_xa__);
disp(sprintf(' %% A_tru_xx__: ')); disp(A_tru_xx__);
disp(sprintf(' %% A_est_xx__: ')); disp(A_est_xx__);
end;%if (flag_verbose>0);
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
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
] = ...
SDE_nlp_ijXaABYC_update_1( ...
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
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_tru_omega ...
,C_tru_l0 ...
,C_tru_l1 ...
);
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% nlp_ijXaABYC_pre %0.6f --> nlp_ijXaABYC_pos %0.6f',parameter_est.nlp_ijXaABYC_pre,parameter_est.nlp_ijXaABYC_pos));
disp(sprintf(' %% nlp_ijXaABYC_integrated_pre %0.6f --> nlp_ijXaABYC_integrated_pos %0.6f',parameter_est.nlp_ijXaABYC_integrated_pre,parameter_est.nlp_ijXaABYC_integrated_pos));
end;%if (flag_verbose>0);
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
xlabel('time');ylabel(tmp_str,'Interpreter','none'); xlim([dolphin_T_ini,dolphin_T_max]);
%s = surfline_0(tmp_W_xt__(1,:),tmp_W_xt__(2,:),t_t_); set(s,'LineWidth',3);
%xlabel(sprintf('%s(1+0,:)',tmp_str),'Interpreter','none');
%ylabel(sprintf('%s(1+1,:)',tmp_str),'Interpreter','none');
title(sprintf('ni %d',ni),'Interpreter','none');
end;%for pcol=0:p_col-1;
end;%for prow=0:p_row-1;
sgtitle(sprintf('X_est_ixt___'),'Interpreter','none');
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
%xlabel('time');ylabel(tmp_str,'Interpreter','none'); xlim([dolphin_T_ini,dolphin_T_max]);
s = surfline_0(tmp_W_xt__(1,:),tmp_W_xt__(2,:),t_t_); set(s,'LineWidth',3);
xlabel(sprintf('%s(1+0,:)',tmp_str),'Interpreter','none');
ylabel(sprintf('%s(1+1,:)',tmp_str),'Interpreter','none');
title(sprintf('ni %d',ni),'Interpreter','none');
end;%for pcol=0:p_col-1;
end;%for prow=0:p_row-1;
sgtitle(sprintf('X_est_ixt___'),'Interpreter','none');
end;%if flag_disp>0;
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
parameter_est.flag_verbose = flag_verbose-2;
parameter_est.flag_disp = flag_verbose-2;
parameter_est.flag_regularize_eccentricity_BtBn = 1;
parameter_est.flag_regularize_eccentricity_simultaneous = 0;
parameter_est.n_iteration_BtBn = 16;
parameter_est.MaxFunEvals_use_BTBn = 1024;
parameter_est.MaxFunEvals_use_simultaneous = dolphin_MaxFunEvals_use_simultaneous;
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
xlabel('time');ylabel(tmp_str,'Interpreter','none'); xlim([dolphin_T_ini,dolphin_T_max]);
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
%xlabel('time');ylabel(tmp_str,'Interpreter','none'); xlim([dolphin_T_ini,dolphin_T_max]);
s = surfline_0(tmp_W_xt__(1,:),tmp_W_xt__(2,:),t_t_); set(s,'LineWidth',3);
xlabel(sprintf('%s(1+0,:)',tmp_str),'Interpreter','none');
ylabel(sprintf('%s(1+1,:)',tmp_str),'Interpreter','none');
title(sprintf('ni %d',ni),'Interpreter','none');
end;%for pcol=0:p_col-1;
end;%for prow=0:p_row-1;
sgtitle(sprintf('all: X_est_ixt___'),'Interpreter','none');
end;%if flag_disp>0;
%%%%%%%%;

if (flag_verbose>0);
disp(sprintf(' %% corr(   A_tru_xx__(:),   A_est_xx__(:)): %+0.6f',corr(A_tru_xx__(:),A_est_xx__(:))));
disp(sprintf(' %% corr(BtBn_tru_xx__(:),BtBn_est_xx__(:)): %+0.6f',corr(BtBn_tru_xx__(:),BtBn_est_xx__(:))));
disp(sprintf(' %% corr(CtCn_tru_xx__(:),CtCn_est_xx__(:)): %+0.6f',corr(CtCn_tru_xx__(:),CtCn_est_xx__(:))));
end;%if (flag_verbose>0); 
