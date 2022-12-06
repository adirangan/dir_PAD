function ...
[ ...
 parameter ...
,nlp_jXaABYC_integrated ...
,X_opt_xt__ ...
,nlp_jXaABYC ...
,nlp_dtXaAB ...
,nlp_jXYC ...
,nlp_dtXaAB_dX_xt__ ...
,nlp_dtXaAB_dX_xtxt__ ...
,nlp_jXYC_dX_xt__ ...
,nlp_jXYC_dX_xtxt__ ...
,BtBn_xx__ ...
,CtCn_xx__ ...
,Q_xt__ ...
,P_xt__ ...
,dX_xdt__ ...
,dQ_xdt__ ...
,dP_xdt__ ...
,ZX_xdt__ ...
,ZQ_xdt__ ...
,ZP_xdt__ ...
] = ...
SDE_nlp_jXaABYC_1( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
,Q_xt__ ...
,P_xt__ ...
,dX_xdt__ ...
,dQ_xdt__ ...
,dP_xdt__ ...
,ZX_xdt__ ...
,ZQ_xdt__ ...
,ZP_xdt__ ...
);

str_thisfunction = 'SDE_nlp_jXaABYC_1';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
nf=0;
rng(0);
n_t = 8; t_t_ = sort(rand(n_t,1));
n_x = 2; X_xt__ = randn(n_x,n_t);
n_j = ceil(n_t*2);
index_nt_from_nj_ = max(0,min(n_t-1,floor(n_t*rand(n_j,1))));
ignore_Y_xj__ = randn(n_x,n_j)>0.5;
Y_xj__ = randn(n_x,n_j);
n_a = 2; a_xa__ = 1*randn(n_x,n_a);
A_xx__ = 2*randn(n_x,n_x);
B_omega = pi/7 ; B_l0 = 1*randn(); B_l1 = 1*randn();
C_omega = pi/5 ; C_l0 = 1*randn(); C_l1 = 1*randn();
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
[ ...
 parameter ...
,nlp_jXaABYC_integrated ...
,X_opt_xt__ ...
,nlp_jXaABYC ...
,nlp_dtXaAB ...
,nlp_jXYC ...
,nlp_dtXaAB_dX_xt__ ...
,nlp_dtXaAB_dX_xtxt__ ...
,nlp_jXYC_dX_xt__ ...
,nlp_jXYC_dX_xtxt__ ...
,BtBn_xx__ ...
,CtCn_xx__ ...
,Q_xt__ ...
,P_xt__ ...
,dX_xdt__ ...
,dQ_xdt__ ...
,dP_xdt__ ...
,ZX_xdt__ ...
,ZQ_xdt__ ...
,ZP_xdt__ ...
] = ...
SDE_nlp_jXaABYC_1( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
disp(sprintf(' %% nlp_jXaABYC_integrated: %f',nlp_jXaABYC_integrated));
disp(sprintf(' %% nlp_jXaABYC: %f',nlp_jXaABYC));
disp(sprintf(' %% nlp_dtXaAB: %f',nlp_dtXaAB));
disp(sprintf(' %% nlp_jXYC: %f',nlp_jXYC));
parameter.flag_verbose=0;
%%%%;
tmp_eps = fnorm(X_xt__)*1e-4;
n_test = 16;
nlp_dtXaAB_dX_emp_l_ = zeros(n_test,1);
nlp_dtXaAB_dX_est_l_ = zeros(n_test,1);
nlp_jXYC_dX_emp_l_ = zeros(n_test,1);
nlp_jXYC_dX_est_l_ = zeros(n_test,1);
for ntest=0:n_test-1;
dX_xt__ = randn(n_x,n_t); dX_xt__ = dX_xt__/fnorm(dX_xt__);
tmp_X_xt__ = X_xt__ + tmp_eps*dX_xt__;
[ ...
 ~ ...
,nlp_jXaABYC_integrated_pos ...
,~ ...
,nlp_jXaABYC_pos ...
,nlp_dtXaAB_pos ...
,nlp_jXYC_pos ...
] = ...
SDE_nlp_jXaABYC_1( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,tmp_X_xt__ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
tmp_X_xt__ = X_xt__ - tmp_eps*dX_xt__;
[ ...
 ~ ...
,nlp_jXaABYC_integrated_neg ...
,~ ...
,nlp_jXaABYC_neg ...
,nlp_dtXaAB_neg ...
,nlp_jXYC_neg ...
] = ...
SDE_nlp_jXaABYC_1( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,tmp_X_xt__ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
nlp_dtXaAB_dX_emp_l_(1+ntest) = (nlp_dtXaAB_pos - nlp_dtXaAB_neg)/(2*tmp_eps);
nlp_dtXaAB_dX_est_l_(1+ntest) = + sum(nlp_dtXaAB_dX_xt__ .* dX_xt__,'all') + reshape(X_xt__,[1,n_x*n_t])*nlp_dtXaAB_dX_xtxt__*reshape(dX_xt__,[n_x*n_t,1]);
nlp_jXYC_dX_emp_l_(1+ntest) = (nlp_jXYC_pos - nlp_jXYC_neg)/(2*tmp_eps);
nlp_jXYC_dX_est_l_(1+ntest) = + sum(nlp_jXYC_dX_xt__ .* dX_xt__,'all') + reshape(X_xt__,[1,n_x*n_t])*nlp_jXYC_dX_xtxt__*reshape(dX_xt__,[n_x*n_t,1]);
end;%for ntest=0:n_test-1;
disp(sprintf(' %% nlp_dtXaAB_dX_emp_l_ vs nlp_dtXaAB_dX_est_l_: %0.16f',fnorm(nlp_dtXaAB_dX_emp_l_ - nlp_dtXaAB_dX_est_l_)/fnorm(nlp_dtXaAB_dX_emp_l_)));
disp(sprintf(' %% nlp_jXYC_dX_emp_l_ vs nlp_jXYC_dX_est_l_: %0.16f',fnorm(nlp_jXYC_dX_emp_l_ - nlp_jXYC_dX_est_l_)/fnorm(nlp_jXYC_dX_emp_l_)));
%%%%;
tmp_RHS_xt_ = reshape(nlp_dtXaAB_dX_xt__,[n_x*n_t,1]) + reshape(nlp_jXYC_dX_xt__,[n_x*n_t,1]) ;
tmp_LHS_xtxt__ = nlp_dtXaAB_dX_xtxt__ + nlp_jXYC_dX_xtxt__ ;
tolerance_master = 1e-6;
%X_opt_xt__ = -reshape(pinv(tmp_LHS_xtxt__,tolerance_master)*tmp_RHS_xt_,[n_x,n_t]);
X_opt_xt__ = -reshape(tmp_LHS_xtxt__\tmp_RHS_xt_,[n_x,n_t]);
tmp_eps = fnorm(X_opt_xt__)*1e-4;
n_test = 16;
scan_max = 0.1;
n_scan = 8+1; eps_s_ = scan_max*transpose(linspace(-1,+1,n_scan));
nlp_jXaABYC_integrated_sl__ = zeros(n_scan,n_test);
nlp_jXaABYC_sl__ = zeros(n_scan,n_test);
for ntest=0:n_test-1;
dX_xt__ = randn(n_x,n_t); dX_xt__ = dX_xt__/fnorm(dX_xt__);
for nscan=0:n_scan-1;
tmp_eps = eps_s_(1+nscan);
tmp_X_xt__ = X_opt_xt__ + tmp_eps*fnorm(X_opt_xt__)*dX_xt__;
[ ...
 ~ ...
,nlp_jXaABYC_integrated_sl__(1+nscan,1+ntest) ...
,~ ...
,nlp_jXaABYC_sl__(1+nscan,1+ntest) ...
] = ...
SDE_nlp_jXaABYC_1( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,tmp_X_xt__ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
end;%for nscan=0:n_scan-1;
end;%for ntest=0:n_test-1;
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,3,1);
plot(nlp_dtXaAB_dX_emp_l_,nlp_dtXaAB_dX_est_l_,'ko');
xlabel('nlp_dtXaAB_dX_emp_l_','Interpreter','none');
ylabel('nlp_dtXaAB_dX_est_l_','Interpreter','none');
title('nlp_dtXaAB','Interpreter','none');
subplot(1,3,2);
plot(nlp_jXYC_dX_emp_l_,nlp_jXYC_dX_est_l_,'ko');
xlabel('nlp_jXYC_dX_emp_l_','Interpreter','none');
ylabel('nlp_jXYC_dX_est_l_','Interpreter','none');
title('nlp_jXYC','Interpreter','none');
subplot(1,3,3);
hold on;
plot(eps_s_,nlp_jXaABYC_sl__,'ko');
plot(eps_s_(1+(n_scan-1)/2),nlp_jXaABYC_sl__(1+(n_scan-1)/2,:),'ko','MarkerFaceColor','r');
hold off;
xlim(scan_max*[-1,+1]); xlabel('scan');
ylabel('nlp_jXaABYC','Interpreter','none');
title('nlp_jXaABYC','Interpreter','none');
%%%%;
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_t=[]; end; na=na+1;
if (nargin<1+na); t_t_=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); X_xt__=[]; end; na=na+1;
if (nargin<1+na); n_a=[]; end; na=na+1;
if (nargin<1+na); a_xa__=[]; end; na=na+1;
if (nargin<1+na); A_xx__=[]; end; na=na+1;
if (nargin<1+na); B_omega=[]; end; na=na+1;
if (nargin<1+na); B_l0=[]; end; na=na+1;
if (nargin<1+na); B_l1=[]; end; na=na+1;
if (nargin<1+na); n_j_=[]; end; na=na+1;
if (nargin<1+na); index_nt_from_nj_=[]; end; na=na+1;
if (nargin<1+na); ignore_Y_xj__=[]; end; na=na+1;
if (nargin<1+na); Y_xj__=[]; end; na=na+1;
if (nargin<1+na); C_omega=[]; end; na=na+1;
if (nargin<1+na); C_l0=[]; end; na=na+1;
if (nargin<1+na); C_l1=[]; end; na=na+1;
if (nargin<1+na); Q_xt__=[]; end; na=na+1;
if (nargin<1+na); P_xt__=[]; end; na=na+1;
if (nargin<1+na); dX_xdt__=[]; end; na=na+1;
if (nargin<1+na); dQ_xdt__=[]; end; na=na+1;
if (nargin<1+na); dP_xdt__=[]; end; na=na+1;
if (nargin<1+na); ZX_xdt__=[]; end; na=na+1;
if (nargin<1+na); ZQ_xdt__=[]; end; na=na+1;
if (nargin<1+na); ZP_xdt__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
if ~isfield(parameter,'tolerance_dt'); parameter.tolerance_dt = 1e-12; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'flag_regularize_eccentricity'); parameter.flag_regularize_eccentricity = 0; end;
tolerance_master = parameter.tolerance_master;
tolerance_dt = parameter.tolerance_dt;
flag_verbose = parameter.flag_verbose;
flag_regularize_eccentricity = parameter.flag_regularize_eccentricity;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%;
tmp_t = tic();
[ ...
 n_nj_from_nt_ ...
,index_nj_from_nt__ ...
] = ...
SDE_index_nj_from_nt_0( ...
 n_t ...
,n_j ...
,index_nt_from_nj_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% SDE_index_nj_from_nt_0: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'SDE_index_nj_from_nt_0',tmp_t,1,max(n_j,n_t));
%%%%;

%%%%;
tmp_t = tic();
[ ...
 nlp_jXaABYC ...
,nlp_dtXaAB ...
,nlp_jXYC ...
,BtBn_xx__ ...
,CtCn_xx__ ...
,Q_xt__ ...
,P_xt__ ...
,dX_xdt__ ...
,dQ_xdt__ ...
,dP_xdt__ ...
,ZX_xdt__ ...
,ZQ_xdt__ ...
,ZP_xdt__ ...
] = ...
SDE_nlp_jXaABYC_strip_1( ...
 n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
,Q_xt__ ...
,P_xt__ ...
,dX_xdt__ ...
,dQ_xdt__ ...
,dP_xdt__ ...
,ZX_xdt__ ...
,ZQ_xdt__ ...
,ZP_xdt__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% SDE_nlp_jXaABYC_strip_1: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'SDE_nlp_jXaABYC_strip_1',tmp_t,1,n_x*n_a*n_t);
%%%%;
[ ...
 ~ ...
,BtBn_xx__ ...
] = ...
SDE_BtBn_0( ...
 [] ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
[ ...
 ~ ...
,CtCn_xx__ ...
] = ...
SDE_BtBn_0( ...
 [] ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
%%%%;

%%%%%%%%%%%%%%%%;
if (n_t<=1);
nlp_jXaABYC_integrated = 0;
X_opt_xt__ = zeros(n_x,n_t);
for nt=0:n_t-1;
tmp_index_ = efind(index_nt_from_nj_==nt);
if numel(tmp_index_)> 0;
nj = min(tmp_index_);
X_opt_xt__(:,1+nt) = Y_xj__(:,1+nj);
end;%if numel(tmp_index_)> 0;
end;%for nt=0:n_t-1;
nlp_jXaABYC = 0;
nlp_dtXaAB = 0;
nlp_jXYC = 0;
nlp_dtXaAB_dX_xt__ = sparse(n_x,n_t);
nlp_dtXaAB_dX_xtxt__ = sparse(n_x*n_t,n_x*n_t);
nlp_jXYC_dX_xt__ = sparse(n_x,n_t);
nlp_jXYC_dX_xtxt__ = sparse(n_x*n_t,n_x*n_t);
end;%if (n_t<=1);
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
if (n_t> 1);
%%%%%%%%%%%%%%%%;

XY_xj__ = Y_xj__ - X_xt__(:,1+index_nt_from_nj_);

if isempty(ignore_Y_xj__); ignore_Y_xj__ = ~isfinite(Y_xj__); end;%if isempty(ignore_Y_xj__); 
ignore_Y_sum_j_ = sum(ignore_Y_xj__,1);
index_use_01_j_ = efind(ignore_Y_sum_j_==0); n_01 = numel(index_use_01_j_);
index_use_0_j_ = efind( (ignore_Y_sum_j_==1) & (ignore_Y_xj__(1+1,:)==1) ); n_0 = numel(index_use_0_j_); %<-- use 0 and ignore 1. ;
index_use_1_j_ = efind( (ignore_Y_sum_j_==1) & (ignore_Y_xj__(1+0,:)==1) ); n_1 = numel(index_use_1_j_); %<-- use 1 and ignore 0. ;
[Z2_base_0,l2_stretch_0] = SDE_missing_2d_integral_helper_0(0,C_omega,C_l0,C_l1); %<-- ignore XY_(1+0) ;
[Z2_base_1,l2_stretch_1] = SDE_missing_2d_integral_helper_0(1,C_omega,C_l0,C_l1); %<-- ignore XY_(1+1) ;
tmp_d_01 = + 0.5*n_x*log(2*pi) - 0.5*(C_l0 + C_l1); %<-- use both. ;
tmp_d_0  = + 0.5*  1*log(2*pi) - 0.5*(C_l0 + C_l1) + 0.5*log(l2_stretch_1); %<-- ignore 1. ;
tmp_d_1  = + 0.5*  1*log(2*pi) - 0.5*(C_l0 + C_l1) + 0.5*log(l2_stretch_0); %<-- ignore 0. ;

%%%%;
tmp_t = tic();
%%%%;
nlp_jXYC_dX_xt__ = sparse(n_x,n_t);
tmp_n_nz = 2*n_01 + 1*n_0 + 1*n_1;
tmp_index_row_nz_ = zeros(tmp_n_nz,1);
tmp_index_col_nz_ = zeros(tmp_n_nz,1);
tmp_val_nz_ = zeros(tmp_n_nz,1);
tmp_nnz=0;
tmp_CtCn_Y_xj__ = CtCn_xx__*Y_xj__(:,1+index_use_01_j_);
tmp_index_row_nz_(1+tmp_nnz+[0:n_01-1]) = 0;
tmp_index_col_nz_(1+tmp_nnz+[0:n_01-1]) = index_nt_from_nj_(1+index_use_01_j_);
tmp_val_nz_(1+tmp_nnz+[0:n_01-1]) = -tmp_CtCn_Y_xj__(1+0,:);
tmp_nnz = tmp_nnz+n_01;
tmp_index_row_nz_(1+tmp_nnz+[0:n_01-1]) = 1;
tmp_index_col_nz_(1+tmp_nnz+[0:n_01-1]) = index_nt_from_nj_(1+index_use_01_j_);
tmp_val_nz_(1+tmp_nnz+[0:n_01-1]) = -tmp_CtCn_Y_xj__(1+1,:);
tmp_nnz = tmp_nnz+n_01;
tmp_Z2_base_1_Y_j_ = Z2_base_1*Y_xj__(1+0,1+index_use_0_j_);
tmp_index_row_nz_(1+tmp_nnz+[0:n_0-1]) = 0;
tmp_index_col_nz_(1+tmp_nnz+[0:n_0-1]) = index_nt_from_nj_(1+index_use_0_j_);
tmp_val_nz_(1+tmp_nnz+[0:n_0-1]) = -tmp_Z2_base_1_Y_j_;
tmp_nnz = tmp_nnz+n_0;
tmp_Z2_base_0_Y_j_ = Z2_base_0*Y_xj__(1+1,1+index_use_1_j_);
tmp_index_row_nz_(1+tmp_nnz+[0:n_1-1]) = 1;
tmp_index_col_nz_(1+tmp_nnz+[0:n_1-1]) = index_nt_from_nj_(1+index_use_1_j_);
tmp_val_nz_(1+tmp_nnz+[0:n_1-1]) = -tmp_Z2_base_0_Y_j_;
tmp_nnz = tmp_nnz+n_1;
nlp_jXYC_dX_xt__ = sparse(1+tmp_index_row_nz_,1+tmp_index_col_nz_,tmp_val_nz_,n_x,n_t);
%%%%;
tmp_n_nz = 4*n_01 + 1*n_0 + 1*n_1;
tmp_index_row_nz_ = zeros(tmp_n_nz,1);
tmp_index_col_nz_ = zeros(tmp_n_nz,1);
tmp_val_nz_ = zeros(tmp_n_nz,1);
tmp_nnz=0;
tmp_index_row_nz_(1+tmp_nnz+[0:n_01-1]) = 0+index_nt_from_nj_(1+index_use_01_j_)*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_01-1]) = 0+index_nt_from_nj_(1+index_use_01_j_)*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_01-1]) = +CtCn_xx__(1+0,1+0);
tmp_nnz = tmp_nnz+n_01;
tmp_index_row_nz_(1+tmp_nnz+[0:n_01-1]) = 0+index_nt_from_nj_(1+index_use_01_j_)*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_01-1]) = 1+index_nt_from_nj_(1+index_use_01_j_)*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_01-1]) = +CtCn_xx__(1+0,1+1);
tmp_nnz = tmp_nnz+n_01;
tmp_index_row_nz_(1+tmp_nnz+[0:n_01-1]) = 1+index_nt_from_nj_(1+index_use_01_j_)*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_01-1]) = 0+index_nt_from_nj_(1+index_use_01_j_)*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_01-1]) = +CtCn_xx__(1+1,1+0);
tmp_nnz = tmp_nnz+n_01;
tmp_index_row_nz_(1+tmp_nnz+[0:n_01-1]) = 1+index_nt_from_nj_(1+index_use_01_j_)*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_01-1]) = 1+index_nt_from_nj_(1+index_use_01_j_)*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_01-1]) = +CtCn_xx__(1+1,1+1);
tmp_nnz = tmp_nnz+n_01;
tmp_index_row_nz_(1+tmp_nnz+[0:n_0-1]) = 0+index_nt_from_nj_(1+index_use_0_j_)*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_0-1]) = 0+index_nt_from_nj_(1+index_use_0_j_)*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_0-1]) = +Z2_base_1;
tmp_nnz = tmp_nnz + n_0;
tmp_index_row_nz_(1+tmp_nnz+[0:n_1-1]) = 1+index_nt_from_nj_(1+index_use_1_j_)*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_1-1]) = 1+index_nt_from_nj_(1+index_use_1_j_)*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_1-1]) = +Z2_base_0;
tmp_nnz = tmp_nnz + n_1;
nlp_jXYC_dX_xtxt__ = sparse(1+tmp_index_row_nz_,1+tmp_index_col_nz_,tmp_val_nz_,n_x*n_t,n_x*n_t);
%%%%;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% nlp_jXYC_dX_xtxt__: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'nlp_jXYC_dX_xtxt__',tmp_t,1,n_x*n_t);
%%%%;

flag_check=0;
if flag_check;
%%%%%%%%;
% vectorize this later. ;
%%%%%%%%;
nlp_jXYC_dX_bkp_xt__ = nlp_jXYC_dX_xt__;
nlp_jXYC_dX_bkp_xtxt__ = nlp_jXYC_dX_xtxt__;
%%%%;
tmp_t = tic();
nlp_jXYC_dX_xt__ = sparse(n_x,n_t);
nlp_jXYC_dX_xtxt__ = sparse(n_x*n_t,n_x*n_t);
tmp_d = + 0.5*n_x*log(2*pi) - 0.5*(C_l0 + C_l1);
for nt=0:n_t-1;
index_nj_from_nt_ = index_nj_from_nt__{1+nt};
n_nj_from_nt = n_nj_from_nt_(1+nt);
assert(n_nj_from_nt==numel(index_nj_from_nt_));
if (n_nj_from_nt> 0);
X_x_ = X_xt__(:,1+nt);
Y_sub_xj__ = Y_xj__(:,1+index_nj_from_nt_);
pXY_sub_xj__ = bsxfun(@minus,Y_sub_xj__,X_x_);
ignore_Y_sub_xj__ = ignore_Y_xj__(:,1+index_nj_from_nt_);
for nnj_from_nt=0:n_nj_from_nt-1;
Y_sub_x_ = Y_sub_xj__(:,1+nnj_from_nt);
ignore_Y_sub_x_ = ignore_Y_sub_xj__(:,1+nnj_from_nt);
if (sum(ignore_Y_sub_x_)==0);
nlp_jXYC_dX_xt__(:,1+nt) = nlp_jXYC_dX_xt__(:,1+nt) + -CtCn_xx__*Y_sub_x_;
nlp_jXYC_dX_xtxt__(1+0+nt*n_x,1+0+nt*n_x) = nlp_jXYC_dX_xtxt__(1+0+nt*n_x,1+0+nt*n_x) + +CtCn_xx__(1+0,1+0);
nlp_jXYC_dX_xtxt__(1+0+nt*n_x,1+1+nt*n_x) = nlp_jXYC_dX_xtxt__(1+0+nt*n_x,1+1+nt*n_x) + +CtCn_xx__(1+0,1+1);
nlp_jXYC_dX_xtxt__(1+1+nt*n_x,1+0+nt*n_x) = nlp_jXYC_dX_xtxt__(1+1+nt*n_x,1+0+nt*n_x) + +CtCn_xx__(1+1,1+0);
nlp_jXYC_dX_xtxt__(1+1+nt*n_x,1+1+nt*n_x) = nlp_jXYC_dX_xtxt__(1+1+nt*n_x,1+1+nt*n_x) + +CtCn_xx__(1+1,1+1);
else;
if (ignore_Y_sub_x_(1+0)==0) & (ignore_Y_sub_x_(1+1)==1) ; %<-- ignore XY_x_(1+1) ;
nlp_jXYC_dX_xt__(1+0,1+nt) = nlp_jXYC_dX_xt__(1+0,1+nt) + -Z2_base_1*Y_sub_x_(1+0);
nlp_jXYC_dX_xtxt__(1+0+nt*n_x,1+0+nt*n_x) = nlp_jXYC_dX_xtxt__(1+0+nt*n_x,1+0+nt*n_x) + +Z2_base_1;
end;%if (ignore_Y_sub_x_(1+0)==0) & (ignore_Y_sub_x_(1+1)==1) ;
if (ignore_Y_sub_x_(1+0)==1) & (ignore_Y_sub_x_(1+1)==0) ; %<-- ignore XY_x_(1+0) ;
nlp_jXYC_dX_xt__(1+1,1+nt) = nlp_jXYC_dX_xt__(1+1,1+nt) + -Z2_base_0*Y_sub_x_(1+1);
nlp_jXYC_dX_xtxt__(1+1+nt*n_x,1+1+nt*n_x) = nlp_jXYC_dX_xtxt__(1+1+nt*n_x,1+1+nt*n_x) + +Z2_base_0;
end;%if (ignore_Y_sub_x_(1+0)==1) & (ignore_Y_sub_x_(1+1)==0) ;
end;%for nnj_from_nt=0:n_nj_from_nt-1;
end;%else;
end;%if (n_nj_from_nt> 0);
end;%for nt=0:n_t-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% nlp_jXYC_dX_xtxt__ (slow): %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'nlp_jXYC_dX_xtxt__ (slow)',tmp_t,1,n_x*n_t);
%%%%;
disp(sprintf(' %% nlp_jXYC_dX_bkp_xt__ vs nlp_jXYC_dX_xt__: %0.16f',fnorm(nlp_jXYC_dX_bkp_xt__-nlp_jXYC_dX_xt__)/fnorm(nlp_jXYC_dX_bkp_xt__)));
disp(sprintf(' %% nlp_jXYC_dX_bkp_xtxt__ vs nlp_jXYC_dX_xtxt__: %0.16f',fnorm(nlp_jXYC_dX_bkp_xtxt__-nlp_jXYC_dX_xtxt__)/fnorm(nlp_jXYC_dX_bkp_xtxt__)));
%%%%%%%%;
end;%if flag_check;

%%%%;
dt_dt_ = diff(t_t_); n_dt = numel(dt_dt_);
dt_min = min(dt_dt_);
if (dt_min< tolerance_dt);
disp(sprintf(' %% Warning, dt_min %0.16f below tolerance_dt %0.16f',dt_min,tolerance_dt));
end;%if (dt_min< tolerance_dt);
M_xx__ = BtBn_xx__;
L_xx__ = transpose(A_xx__)*BtBn_xx__;
R_xx__ = BtBn_xx__*A_xx__;
G_xx__ = transpose(A_xx__)*BtBn_xx__*A_xx__;
%%%%;

%%%%;
tmp_t = tic();
%%%%;
Q_pos_xdt__ = Q_xt__(:,2:n_t-0);
Q_pre_xdt__ = Q_xt__(:,1:n_t-1);
MdQ_xdt__ = bsxfun(@rdivide,M_xx__*(-dQ_xdt__),reshape(max(1e-12,dt_dt_),[1,n_dt]));
LQ_xdt__ = L_xx__*(-dQ_xdt__);
RQ_xdt__ = R_xx__*(-Q_pre_xdt__);
GQ_xdt__ = bsxfun(@times,G_xx__*(-Q_pre_xdt__),reshape(max(1e-12,dt_dt_),[1,n_dt]));
nlp_dtXaAB_dX_bot_xt__ = zeros(n_x,n_t);
nlp_dtXaAB_dX_bot_xt__(:,1:n_t-1) = - MdQ_xdt__ - LQ_xdt__ + RQ_xdt__ + GQ_xdt__;
nlp_dtXaAB_dX_top_xt__ = zeros(n_x,n_t);
nlp_dtXaAB_dX_top_xt__(:,2:n_t-0) = + MdQ_xdt__ - RQ_xdt__;
%%%%;
tmp_n_nz = n_dt*8;
tmp_index_row_nz_ = zeros(tmp_n_nz,1);
tmp_index_col_nz_ = zeros(tmp_n_nz,1);
tmp_val_nz_ = zeros(tmp_n_nz,1);
tmp_nnz=0;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[0:n_dt-1]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[1:n_dt-0]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = - M_xx__(1+0,1+0)./max(1e-12,dt_dt_) - L_xx__(1+0,1+0);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[0:n_dt-1]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[1:n_dt-0]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = - M_xx__(1+0,1+1)./max(1e-12,dt_dt_) - L_xx__(1+0,1+1);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[0:n_dt-1]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[1:n_dt-0]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = - M_xx__(1+1,1+0)./max(1e-12,dt_dt_) - L_xx__(1+1,1+0);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[0:n_dt-1]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[1:n_dt-0]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = - M_xx__(1+1,1+1)./max(1e-12,dt_dt_) - L_xx__(1+1,1+1);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[0:n_dt-1]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[0:n_dt-1]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = + M_xx__(1+0,1+0)./max(1e-12,dt_dt_) + L_xx__(1+0,1+0) + R_xx__(1+0,1+0) + G_xx__(1+0,1+0).*max(1e-12,dt_dt_);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[0:n_dt-1]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[0:n_dt-1]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = + M_xx__(1+0,1+1)./max(1e-12,dt_dt_) + L_xx__(1+0,1+1) + R_xx__(1+0,1+1) + G_xx__(1+0,1+1).*max(1e-12,dt_dt_);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[0:n_dt-1]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[0:n_dt-1]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = + M_xx__(1+1,1+0)./max(1e-12,dt_dt_) + L_xx__(1+1,1+0) + R_xx__(1+1,1+0) + G_xx__(1+1,1+0).*max(1e-12,dt_dt_);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[0:n_dt-1]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[0:n_dt-1]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = + M_xx__(1+1,1+1)./max(1e-12,dt_dt_) + L_xx__(1+1,1+1) + R_xx__(1+1,1+1) + G_xx__(1+1,1+1).*max(1e-12,dt_dt_);
tmp_nnz = tmp_nnz + n_dt;
nlp_dtXaAB_dX_bot_xtxt__ = sparse(1+tmp_index_row_nz_,1+tmp_index_col_nz_,tmp_val_nz_,n_x*n_t,n_x*n_t);
%%%%;
tmp_n_nz = n_dt*8;
tmp_index_row_nz_ = zeros(tmp_n_nz,1);
tmp_index_col_nz_ = zeros(tmp_n_nz,1);
tmp_val_nz_ = zeros(tmp_n_nz,1);
tmp_nnz=0;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[1:n_dt-0]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[0:n_dt-1]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = - M_xx__(1+0,1+0)./max(1e-12,dt_dt_) - R_xx__(1+0,1+0);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[1:n_dt-0]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[0:n_dt-1]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = - M_xx__(1+0,1+1)./max(1e-12,dt_dt_) - R_xx__(1+0,1+1);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[1:n_dt-0]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[0:n_dt-1]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = - M_xx__(1+1,1+0)./max(1e-12,dt_dt_) - R_xx__(1+1,1+0);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[1:n_dt-0]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[0:n_dt-1]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = - M_xx__(1+1,1+1)./max(1e-12,dt_dt_) - R_xx__(1+1,1+1);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[1:n_dt-0]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[1:n_dt-0]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = + M_xx__(1+0,1+0)./max(1e-12,dt_dt_);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[1:n_dt-0]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[1:n_dt-0]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = + M_xx__(1+0,1+1)./max(1e-12,dt_dt_);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[1:n_dt-0]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 0+[1:n_dt-0]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = + M_xx__(1+1,1+0)./max(1e-12,dt_dt_);
tmp_nnz = tmp_nnz + n_dt;
tmp_index_row_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[1:n_dt-0]*n_x;
tmp_index_col_nz_(1+tmp_nnz+[0:n_dt-1]) = 1+[1:n_dt-0]*n_x;
tmp_val_nz_(1+tmp_nnz+[0:n_dt-1]) = + M_xx__(1+1,1+1)./max(1e-12,dt_dt_);
tmp_nnz = tmp_nnz + n_dt;
nlp_dtXaAB_dX_top_xtxt__ = sparse(1+tmp_index_row_nz_,1+tmp_index_col_nz_,tmp_val_nz_,n_x*n_t,n_x*n_t);
%%%%;
%%%%;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% nlp_dtXaAB_dX_xtxt__: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'nlp_dtXaAB_dX_xtxt__',tmp_t,1,n_x*n_t);
%%%%;

flag_check=0;
if flag_check;
%%%%%%%%;
% vectorize this later. ;
%%%%%%%%;
nlp_dtXaAB_dX_bot_bkp_xt__ = nlp_dtXaAB_dX_bot_xt__;
nlp_dtXaAB_dX_top_bkp_xt__ = nlp_dtXaAB_dX_top_xt__;
nlp_dtXaAB_dX_bot_bkp_xtxt__ = nlp_dtXaAB_dX_bot_xtxt__;
nlp_dtXaAB_dX_top_bkp_xtxt__ = nlp_dtXaAB_dX_top_xtxt__;
tmp_t = tic();
nlp_dtXaAB_dX_bot_xt__ = sparse(n_x,n_t);
nlp_dtXaAB_dX_bot_xtxt__ = sparse(n_x*n_t,n_x*n_t);
nlp_dtXaAB_dX_top_xt__ = sparse(n_x,n_t);
nlp_dtXaAB_dX_top_xtxt__ = sparse(n_x*n_t,n_x*n_t);
for nt=0:n_t-1;
if (nt< n_t-1);
Q_pos_x_ = Q_xt__(:,1+nt+1);
Q_pre_x_ = Q_xt__(:,1+nt+0);
ndt = nt;
dQ_x_ = dQ_xdt__(:,1+ndt);
tmp_dt = max(1e-12,dt_dt_(1+ndt));
%nlp_dtXaAB_dX_bot_xt__(:,1+nt) = - M_xx__/tmp_dt * (-Q_pos_x_ + Q_pre_x_) - L_xx__*(-Q_pos_x_ + Q_pre_x_) + R_xx__*(-Q_pre_x_) + tmp_dt*G_xx__*(-Q_pre_x_) ;
nlp_dtXaAB_dX_bot_xt__(:,1+nt) = - M_xx__/tmp_dt * (-dQ_x_) - L_xx__*(-dQ_x_) + R_xx__*(-Q_pre_x_) + tmp_dt*G_xx__*(-Q_pre_x_) ;
nt_out = nt; nt_0in = nt_out + 1;
nlp_dtXaAB_dX_bot_xtxt__(1+0+nt_out*n_x,1+0+nt_0in*n_x) = - M_xx__(1+0,1+0)/tmp_dt - L_xx__(1+0,1+0) ;
nlp_dtXaAB_dX_bot_xtxt__(1+0+nt_out*n_x,1+1+nt_0in*n_x) = - M_xx__(1+0,1+1)/tmp_dt - L_xx__(1+0,1+1) ;
nlp_dtXaAB_dX_bot_xtxt__(1+1+nt_out*n_x,1+0+nt_0in*n_x) = - M_xx__(1+1,1+0)/tmp_dt - L_xx__(1+1,1+0) ;
nlp_dtXaAB_dX_bot_xtxt__(1+1+nt_out*n_x,1+1+nt_0in*n_x) = - M_xx__(1+1,1+1)/tmp_dt - L_xx__(1+1,1+1) ;
nt_out = nt; nt_0in = nt_out + 0;
nlp_dtXaAB_dX_bot_xtxt__(1+0+nt_out*n_x,1+0+nt_0in*n_x) = + M_xx__(1+0,1+0)/tmp_dt + L_xx__(1+0,1+0) + R_xx__(1+0,1+0) + tmp_dt*G_xx__(1+0,1+0) ;
nlp_dtXaAB_dX_bot_xtxt__(1+0+nt_out*n_x,1+1+nt_0in*n_x) = + M_xx__(1+0,1+1)/tmp_dt + L_xx__(1+0,1+1) + R_xx__(1+0,1+1) + tmp_dt*G_xx__(1+0,1+1) ;
nlp_dtXaAB_dX_bot_xtxt__(1+1+nt_out*n_x,1+0+nt_0in*n_x) = + M_xx__(1+1,1+0)/tmp_dt + L_xx__(1+1,1+0) + R_xx__(1+1,1+0) + tmp_dt*G_xx__(1+1,1+0) ;
nlp_dtXaAB_dX_bot_xtxt__(1+1+nt_out*n_x,1+1+nt_0in*n_x) = + M_xx__(1+1,1+1)/tmp_dt + L_xx__(1+1,1+1) + R_xx__(1+1,1+1) + tmp_dt*G_xx__(1+1,1+1) ;
end;%if (nt< n_t-1);
if (nt> 0);
Q_pos_x_ = Q_xt__(:,1+nt+0);
Q_pre_x_ = Q_xt__(:,1+nt-1);
ndt = nt-1;
dQ_x_ = dQ_xdt__(:,1+ndt);
tmp_dt = max(1e-12,dt_dt_(1+ndt));
%nlp_dtXaAB_dX_top_xt__(:,1+nt) = + M_xx__/tmp_dt * (-Q_pos_x_ + Q_pre_x_) - R_xx__*(-Q_pre_x_) ;
nlp_dtXaAB_dX_top_xt__(:,1+nt) = + M_xx__/tmp_dt * (-dQ_x_) - R_xx__*(-Q_pre_x_) ;
nt_out = nt; nt_0in = nt_out - 1;
nlp_dtXaAB_dX_top_xtxt__(1+0+nt_out*n_x,1+0+nt_0in*n_x) = - M_xx__(1+0,1+0)/tmp_dt - R_xx__(1+0,1+0) ;
nlp_dtXaAB_dX_top_xtxt__(1+0+nt_out*n_x,1+1+nt_0in*n_x) = - M_xx__(1+0,1+1)/tmp_dt - R_xx__(1+0,1+1) ;
nlp_dtXaAB_dX_top_xtxt__(1+1+nt_out*n_x,1+0+nt_0in*n_x) = - M_xx__(1+1,1+0)/tmp_dt - R_xx__(1+1,1+0) ;
nlp_dtXaAB_dX_top_xtxt__(1+1+nt_out*n_x,1+1+nt_0in*n_x) = - M_xx__(1+1,1+1)/tmp_dt - R_xx__(1+1,1+1) ;
nt_out = nt; nt_0in = nt_out - 0;
nlp_dtXaAB_dX_top_xtxt__(1+0+nt_out*n_x,1+0+nt_0in*n_x) = + M_xx__(1+0,1+0)/tmp_dt ;
nlp_dtXaAB_dX_top_xtxt__(1+0+nt_out*n_x,1+1+nt_0in*n_x) = + M_xx__(1+0,1+1)/tmp_dt ;
nlp_dtXaAB_dX_top_xtxt__(1+1+nt_out*n_x,1+0+nt_0in*n_x) = + M_xx__(1+1,1+0)/tmp_dt ;
nlp_dtXaAB_dX_top_xtxt__(1+1+nt_out*n_x,1+1+nt_0in*n_x) = + M_xx__(1+1,1+1)/tmp_dt ;
end;%if (nt> 0);
end;%for nt=0:n_t-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% nlp_dtXaAB_dX_xtxt__ (slow): %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'nlp_dtXaAB_dX_xtxt__ (slow)',tmp_t,1,n_x*n_t);
%%%%;
disp(sprintf(' %% nlp_dtXaAB_dX_bot_bkp_xt__ vs nlp_dtXaAB_dX_bot_xt__: %0.16f',fnorm(nlp_dtXaAB_dX_bot_bkp_xt__-nlp_dtXaAB_dX_bot_xt__)/fnorm(nlp_dtXaAB_dX_bot_bkp_xt__)));
disp(sprintf(' %% nlp_dtXaAB_dX_top_bkp_xt__ vs nlp_dtXaAB_dX_top_xt__: %0.16f',fnorm(nlp_dtXaAB_dX_top_bkp_xt__-nlp_dtXaAB_dX_top_xt__)/fnorm(nlp_dtXaAB_dX_top_bkp_xt__)));
disp(sprintf(' %% nlp_dtXaAB_dX_bot_bkp_xtxt__ vs nlp_dtXaAB_dX_bot_xtxt__: %0.16f',fnorm(nlp_dtXaAB_dX_bot_bkp_xtxt__-nlp_dtXaAB_dX_bot_xtxt__)/fnorm(nlp_dtXaAB_dX_bot_bkp_xtxt__)));
disp(sprintf(' %% nlp_dtXaAB_dX_top_bkp_xtxt__ vs nlp_dtXaAB_dX_top_xtxt__: %0.16f',fnorm(nlp_dtXaAB_dX_top_bkp_xtxt__-nlp_dtXaAB_dX_top_xtxt__)/fnorm(nlp_dtXaAB_dX_top_bkp_xtxt__)));
%%%%%%%%;
end;%if flag_check;

nlp_dtXaAB_dX_xt__ = nlp_dtXaAB_dX_bot_xt__ + nlp_dtXaAB_dX_top_xt__;
nlp_dtXaAB_dX_xtxt__ = nlp_dtXaAB_dX_bot_xtxt__ + nlp_dtXaAB_dX_top_xtxt__;

%%%%;
tmp_RHS_xt_ = reshape(nlp_dtXaAB_dX_xt__,[n_x*n_t,1]) + reshape(nlp_jXYC_dX_xt__,[n_x*n_t,1]) ;
tmp_LHS_xtxt__ = nlp_dtXaAB_dX_xtxt__ + nlp_jXYC_dX_xtxt__ ;
if n_t< 128;
X_opt_xt__ = -reshape(pinv(full(tmp_LHS_xtxt__),tolerance_master)*tmp_RHS_xt_,[n_x,n_t]);
end;%if n_t< 128;
if n_t>=128;
X_opt_xt__ = -reshape(tmp_LHS_xtxt__\tmp_RHS_xt_,[n_x,n_t]);
end;%if n_t>=128;
tmp_t = tic();
[tmp_LHS_l_xtxt__,tmp_LHS_u_xtxt__,tmp_LHS_t_xtxt__] = lu(tmp_LHS_xtxt__);
tmp_LHS_ldetu = sum(log(abs(diag(tmp_LHS_u_xtxt__))));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% tmp_LHS_ldetu: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'tmp_LHS_ldetu',tmp_t,1,nnz(tmp_LHS_xtxt__));
nlp_jXaABYC_integrated = nlp_dtXaAB + nlp_jXYC - 0.5*n_t*n_x*log(2*pi) + 0.5*tmp_LHS_ldetu; %<-- multidimensional laplace integral. ;
%%%%;

%%%%%%%%%%%%%%%%;
end;%if (n_t> 1);
%%%%%%%%%%%%%%%%;

nlp_jXaABYC_integrated = nlp_jXaABYC_integrated ...
+ flag_regularize_eccentricity* ...
SDE_nlp_jXYC_eccentricity_1( ...
 n_x ...
,n_dt ...
,B_l0 ...
,B_l1 ...
);
nlp_jXaABYC_integrated = nlp_jXaABYC_integrated ...
+ flag_regularize_eccentricity* ...
SDE_nlp_jXYC_eccentricity_1( ...
 n_x ...
,n_j ...
,C_l0 ...
,C_l1 ...
);


if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;





