function ...
[ ...
 parameter ...
,nlp ...
,nlp_tXaAB ...
,nlp_tXYC ...
,nlp_tXYC_dX_xt__ ...
,nlp_tXYC_dX_xtxt____ ...
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
PAD_nlp_tXaABYC_0( ...
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
,ignore_Y_xt__ ...
,Y_xt__ ...
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

str_thisfunction = 'PAD_nlp_tXaABYC_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
nf=0;
rng(0);
n_t = 8; t_t_ = sort(rand(n_t,1));
n_x = 2; X_xt__ = randn(n_x,n_t);
Y_xt__ = randn(n_x,n_t); ignore_Y_xt__ = rand(n_x,n_t)< 0.5;
n_a = 2; a_xa__ = randn(n_x,n_a);
A_xx__ = randn(n_x,n_x);
B_omega = pi/7 ; B_l0 = randn(); B_l1 = randn();
C_omega = pi/5 ; C_l0 = randn(); C_l1 = randn();
parameter = struct('type','parameter');
[ ...
 parameter ...
,nlp ...
,nlp_tXaAB ...
,nlp_tXYC ...
,nlp_tXYC_dX_xt__ ...
,nlp_tXYC_dX_xtxt____ ...
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
PAD_nlp_tXaABYC_0( ...
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
,ignore_Y_xt__ ...
,Y_xt__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
disp(sprintf(' %% nlp: %f',nlp));
disp(sprintf(' %% nlp_tXaAB: %f',nlp_tXaAB));
disp(sprintf(' %% nlp_tXYC: %f',nlp_tXYC));
%%%%;
tmp_eps = fnorm(X_xt__)*1e-4;
n_test = 16;
nlp_tXYC_dX_emp_l_ = zeros(n_test,1);
nlp_tXYC_dX_est_l_ = zeros(n_test,1);
for ntest=0:n_test-1;
dX_xt__ = randn(n_x,n_t); dX_xt__ = dX_xt__/fnorm(dX_xt__);
tmp_X_xt__ = X_xt__ + tmp_eps*dX_xt__;
[ ...
 ~ ...
,nlp_pos ...
,nlp_tXaAB_pos ...
,nlp_tXYC_pos ...
] = ...
PAD_nlp_tXaABYC_0( ...
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
,ignore_Y_xt__ ...
,Y_xt__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
tmp_X_xt__ = X_xt__ - tmp_eps*dX_xt__;
[ ...
 ~ ...
,nlp_neg ...
,nlp_tXaAB_neg ...
,nlp_tXYC_neg ...
] = ...
PAD_nlp_tXaABYC_0( ...
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
,ignore_Y_xt__ ...
,Y_xt__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
nlp_tXYC_dX_emp_l_(1+ntest) = (nlp_tXYC_pos - nlp_tXYC_neg)/(2*tmp_eps);
nlp_tXYC_dX_est_l_(1+ntest) = sum(nlp_tXYC_dX_xt__ .* dX_xt__,'all') + reshape(X_xt__,[1,n_x*n_t])*nlp_tXYC_dX_xtxt____*reshape(dX_xt__,[n_x*n_t,1]);
end;%for ntest=0:n_test-1;
disp(sprintf(' %% nlp_tXYC_dX_emp_l_ vs nlp_tXYC_dX_est_l_: %0.16f',fnorm(nlp_tXYC_dX_emp_l_ - nlp_tXYC_dX_est_l_)/fnorm(nlp_tXYC_dX_emp_l_)));
%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
plot(nlp_tXYC_dX_emp_l_,nlp_tXYC_dX_est_l_,'ko');
xlabel('nlp_tXYC_dX_emp_l_','Interpreter','none');
ylabel('nlp_tXYC_dX_est_l_','Interpreter','none');
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
if (nargin<1+na); ignore_Y_xt__=[]; end; na=na+1;
if (nargin<1+na); Y_xt__=[]; end; na=na+1;
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
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;

%%%%;
[ ...
 nlp ...
,nlp_tXaAB ...
,nlp_tXYC ...
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
PAD_nlp_tXaABYC_strip_0( ...
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
,ignore_Y_xt__ ...
,Y_xt__ ...
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
%%%%;
[ ...
 ~ ...
,CtCn_xx__ ...
] = ...
PAD_BtBn_0( ...
 [] ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
[Z2_base_0,l2_stretch_0] = PAD_missing_2d_integral_helper_0(0,C_omega,C_l0,C_l1); %<-- ignore XY_(1+0) ;
[Z2_base_1,l2_stretch_1] = PAD_missing_2d_integral_helper_0(1,C_omega,C_l0,C_l1); %<-- ignore XY_(1+1) ;
%%%%;
if isempty(ignore_Y_xt__); ignore_Y_xt__ = ~isfinite(Y_xt__); end;%if isempty(ignore_Y_xt__); 
%%%%;
XY_xt__ = Y_xt__ - X_xt__;
%%%%%%%%;
% vectorize this later. ;
%%%%%%%%;
nlp_tXYC_dX_xt__ = sparse(n_x,n_t);
nlp_tXYC_dX_xtxt____ = sparse(n_x*n_t,n_x*n_t);
tmp_d = + 0.5*n_x*log(2*pi) - 0.5*(C_l0 + C_l1);
for nt=0:n_t-1;
X_x_ = X_xt__(:,1+nt);
Y_x_ = Y_xt__(:,1+nt);
XY_x_ = XY_xt__(:,1+nt);
ignore_Y_x_ = ignore_Y_xt__(:,1+nt);
if (sum(ignore_Y_x_)==0);
nlp_tXYC_dX_xt__(:,1+nt) = -CtCn_xx__*Y_x_;
nlp_tXYC_dX_xtxt____(1+0+nt*n_x,1+0+nt*n_x) = +CtCn_xx__(1+0,1+0);
nlp_tXYC_dX_xtxt____(1+0+nt*n_x,1+1+nt*n_x) = +CtCn_xx__(1+0,1+1);
nlp_tXYC_dX_xtxt____(1+1+nt*n_x,1+0+nt*n_x) = +CtCn_xx__(1+1,1+0);
nlp_tXYC_dX_xtxt____(1+1+nt*n_x,1+1+nt*n_x) = +CtCn_xx__(1+1,1+1);
else;
if (ignore_Y_x_(1+0)==0) & (ignore_Y_x_(1+1)==1) ; %<-- ignore XY_x_(1+1) ;
nlp_tXYC_dX_xt__(1+0,1+nt) = -Z2_base_1*Y_x_(1+0);
nlp_tXYC_dX_xtxt____(1+0+nt*n_x,1+0+nt*n_x) = +Z2_base_1;
end;%if (ignore_Y_x_(1+0)==0) & (ignore_Y_x_(1+1)==1) ;
if (ignore_Y_x_(1+0)==1) & (ignore_Y_x_(1+1)==0) ; %<-- ignore XY_x_(1+0) ;
nlp_tXYC_dX_xt__(1+1,1+nt) = -Z2_base_0*Y_x_(1+1);
nlp_tXYC_dX_xtxt____(1+1+nt*n_x,1+1+nt*n_x) = +Z2_base_0;
end;%if (ignore_Y_x_(1+0)==1) & (ignore_Y_x_(1+1)==0) ;
end;%else;
end;%for nt=0:n_dt-1;
%%%%;






