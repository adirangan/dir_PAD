function ...
[ ...
 parameter ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_ixt___ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,ignore_Y_ixt___ ...
,Y_ixt___ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
] = ...
PAD_nlp_itXaABYC_update_0( ...
 parameter ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_ixt___ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,ignore_Y_ixt___ ...
,Y_ixt___ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);

str_thisfunction = 'PAD_nlp_itXaABYC_update_0';

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
if (nargin<1+na); X_ixt___=[]; end; na=na+1;
if (nargin<1+na); n_a=[]; end; na=na+1;
if (nargin<1+na); a_xa__=[]; end; na=na+1;
if (nargin<1+na); A_xx__=[]; end; na=na+1;
if (nargin<1+na); B_omega=[]; end; na=na+1;
if (nargin<1+na); B_l0=[]; end; na=na+1;
if (nargin<1+na); B_l1=[]; end; na=na+1;
if (nargin<1+na); ignore_Y_ixt___=[]; end; na=na+1;
if (nargin<1+na); Y_ixt___=[]; end; na=na+1;
if (nargin<1+na); C_omega=[]; end; na=na+1;
if (nargin<1+na); C_l0=[]; end; na=na+1;
if (nargin<1+na); C_l1=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'str_update'); parameter.str_update = 'CtCn_xx__'; end;
if ~isfield(parameter,'flag_regularize_eccentricity'); parameter.flag_regularize_eccentricity = 0; end;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;
str_update = parameter.str_update;
flag_regularize_eccentricity = parameter.flag_regularize_eccentricity;

if flag_verbose; disp(sprintf(' %% str_update: %s',str_update)); end;

if isempty(X_ixt___);
X_ixt___ = cell(n_i,1);
for ni=0:n_i-1;
ignore_Y_xt__ = ignore_Y_ixt___{1+ni};
Y_xt__ = Y_ixt___{1+ni};
X_xt__ = 0.5*Y_xt__.*(~ignore_Y_xt__);
X_ixt___{1+ni} = X_xt__;
clear ignore_Y_xt__ Y_xt__ X_xt__ ;
end;%for ni=0:n_i-1;
end;%if isempty(X_ixt___);

if isempty(a_xa__); a_xa__ = zeros(n_x,n_a); end;
if isempty(A_xx__); A_xx__ = zeros(n_x,n_x); end;
if isempty(B_omega); B_omega = 0; end; if isempty(B_l0); B_l0 = 0; end; if isempty(B_l1); B_l1 = 0; end;
if isempty(C_omega); C_omega = 0; end; if isempty(C_l0); C_l0 = 0; end; if isempty(C_l1); C_l1 = 0; end;

if ~isempty(strfind(str_update,'CtCn_xx__'));
%%%%%%%%;
n_t_all = sum(n_t_i_);
X_all_xt__ = zeros(n_x,n_t_all);
ignore_Y_all_xt__ = zeros(n_x,n_t_all);
Y_all_xt__ = zeros(n_x,n_t_all);
nt_sum=0;
for ni=0:n_i-1;
tmp_index_ = nt_sum+[0:n_t_i_(1+ni)-1];
X_all_xt__(:,1+tmp_index_) = X_ixt___{1+ni};
ignore_Y_all_xt__(:,1+tmp_index_) = ignore_Y_ixt___{1+ni};
Y_all_xt__(:,1+tmp_index_) = Y_ixt___{1+ni};
nt_sum=nt_sum + n_t_i_(1+ni);
end;%for ni=0:n_i-1;
assert(nt_sum==n_t_all);
f_nlp = @(w0l0l1_) PAD_nlp_tXYC_strip_1(n_x,n_t_all,X_all_xt__,ignore_Y_all_xt__,Y_all_xt__,w0l0l1_(1+0),w0l0l1_(1+1),w0l0l1_(1+2)) + flag_regularize_eccentricity*PAD_nlp_tXYC_eccentricity_0(n_x,n_t_all,w0l0l1_(1+1),w0l0l1_(1+2));
nlp_tXYC_pre = f_nlp([C_omega,C_l0,C_l1]);
[w0l0l1_opt_,fval_opt] = fminsearch(f_nlp,[C_omega,C_l0,C_l1],optimset('MaxFunEvals',1024));
C_omega = w0l0l1_opt_(1+0); C_l0 = w0l0l1_opt_(1+1); C_l1 = w0l0l1_opt_(1+2);
nlp_tXYC_pos = f_nlp([C_omega,C_l0,C_l1]);
if flag_verbose;
disp(sprintf(' %% C_omega %+0.2f C_l0 %+0.2f C_l1 %+0.2f: nlp_tXYC %0.6f --> %0.6f',C_omega,C_l0,C_l1,nlp_tXYC_pre,nlp_tXYC_pos));
end;%if flag_verbose;
parameter.nlp_tXYC_pre = nlp_tXYC_pre;
parameter.nlp_tXYC_pos = nlp_tXYC_pos;
clear tmp_index_ X_all_xt__ ignore_Y_all_xt__ Y_all_xt__ f_nlp ;
%%%%%%%%;
end;%if ~isempty(strfind(str_update,'CtCn_xx__'));

if ~isempty(strfind(str_update,'BtBn_xx__'));
%%%%%%%%;
% Designed with the assumption that all the dt are the same. ;
%%%%%%%%;
n_dt_all = sum(n_t_i_-1);
XY_all_xdt__ = zeros(n_x,n_dt_all);
dt_all_dt_ = zeros(n_dt_all,1);
ndt=0;
for ni=0:n_i-1;
[ ...
 nlp_dtXaAB ...
,BtBn_xx__ ...
,Q_xt__ ...
,P_xt__ ...
,dX_xdt__ ...
,dQ_xdt__ ...
,dP_xdt__ ...
,ZX_xdt__ ...
,ZQ_xdt__ ...
,ZP_xdt__ ...
] = ...
PAD_nlp_dtXaAB_strip_0( ...
 n_t_i_(1+ni) ...
,t_it__{1+ni} ...
,n_x ...
,X_ixt___{1+ni} ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
ZP_all_xdt__(:,1+ndt+[0:n_t_i_(1+ni)-1-1]) = ZP_xdt__;
dt_all_dt_(1+ndt+[0:n_t_i_(1+ni)-1-1]) = diff(t_it__{1+ni});
ndt=ndt + n_t_i_(1+ni)-1;
clear nlp_dtXaAB BtBn_xx__ Q_xt__ P_xt__ dX_xdt__ dQ_xdt__ dP_xdt__ ZX_xdt__ ZQ_xdt__ ZP_xdt__ ;
end;%for ni=0:n_i-1;
assert(ndt==n_dt_all);
f_nlp = @(w0l0l1_) PAD_nlp_tXYC_strip_0(n_x,n_dt_all,ZP_all_xdt__,w0l0l1_(1+0),w0l0l1_(1+1),w0l0l1_(1+2)) + flag_regularize_eccentricity*PAD_nlp_tXYC_eccentricity_0(n_x,n_dt_all,w0l0l1_(1+1),w0l0l1_(1+2)); ;
nlp_dtZPB_pre = f_nlp([B_omega,B_l0,B_l1]);
[w0l0l1_opt_,fval_opt] = fminsearch(f_nlp,[B_omega,B_l0,B_l1],optimset('MaxFunEvals',1024));
B_omega = w0l0l1_opt_(1+0); B_l0 = w0l0l1_opt_(1+1); B_l1 = w0l0l1_opt_(1+2);
B_l0 = B_l0 + log(mean(dt_all_dt_));
B_l1 = B_l1 + log(mean(dt_all_dt_));
nlp_dtZPB_pos = f_nlp([B_omega,B_l0,B_l1]);
if flag_verbose;
disp(sprintf(' %% B_omega %+0.2f B_l0 %+0.2f B_l1 %+0.2f: nlp_dtZPB %0.6f --> %0.6f',B_omega,B_l0,B_l1,nlp_dtZPB_pre,nlp_dtZPB_pos));
end;%if flag_verbose;
parameter.nlp_dtZPB_pre = nlp_dtZPB_pre;
parameter.nlp_dtZPB_pos = nlp_dtZPB_pos;
clear ZP_all_xdt__ f_nlp ;
%%%%%%%%;
end;%if ~isempty(strfind(str_update,'BtBn_xx__'));

if ~isempty(strfind(str_update,'a_xa_')) | ~isempty(strfind(str_update,'A_xx__'));
%%%%%%%%;
% Note that we actually update a_xa__ and A_xx__ independently ;
% (each is updated assuming that the other is a constant). ;
% In future versions we should update these simultaneously. ;
%%%%%%%%;
nlp_dtXaAB_sum_pre = 0;
nlp_dtXaAB_da_sum_xa__ = zeros(n_x,n_a);
nlp_dtXaAB_da_sum_xaxa____ = zeros(n_x,n_a,n_x,n_a);
nlp_dtXaAB_dA_sum_xx__ = zeros(n_x,n_x);
nlp_dtXaAB_dA_sum_xxxx____ = zeros(n_x,n_x,n_x,n_x);
for ni=0:n_i-1;
n_t = n_t_i_(1+ni);
t_t_ = t_it__{1+ni};
X_xt__ = X_ixt___{1+ni};
[ ...
 parameter ...
,nlp_dtXaAB ...
,BtBn_xx__ ...
,nlp_dtXaAB_0 ...
,nlp_dtXaAB_a ...
,nlp_dtXaAB_da_xa__ ...
,nlp_dtXaAB_da_xaxa____ ...
,nlp_dtXaAB_A ...
,nlp_dtXaAB_dA_xx__ ...
,nlp_dtXaAB_dA_xxxx____ ...
] = ...
PAD_nlp_dtXaAB_0( ...
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
);
nlp_dtXaAB_sum_pre = nlp_dtXaAB_sum_pre + nlp_dtXaAB;
nlp_dtXaAB_da_sum_xa__ = nlp_dtXaAB_da_sum_xa__ + nlp_dtXaAB_da_xa__;
nlp_dtXaAB_da_sum_xaxa____ = nlp_dtXaAB_da_sum_xaxa____ + nlp_dtXaAB_da_xaxa____;
nlp_dtXaAB_dA_sum_xx__ = nlp_dtXaAB_dA_sum_xx__ + nlp_dtXaAB_dA_xx__;
nlp_dtXaAB_dA_sum_xxxx____ = nlp_dtXaAB_dA_sum_xxxx____ + nlp_dtXaAB_dA_xxxx____;
clear n_t t_t_ X_xt_ nlp_dtXaAB BtBn_xx__ nlp_dtXaAB_0 nlp_dtXaAB_a nlp_dtXaAB_da_xa__ nlp_dtXaAB_da_xaxa____ nlp_dtXaAB_A nlp_dtXaAB_dA_xx__ nlp_dtXaAB_dA_xxxx____ ;
end;%for ni=0:n_i-1;
a_xa__ = reshape( - pinv(reshape(nlp_dtXaAB_da_sum_xaxa____,[n_x*n_a,n_x*n_a]),tolerance_master) * reshape(nlp_dtXaAB_da_sum_xa__,[n_x*n_a,1]) , [n_x,n_a] );
A_xx__ = reshape( - pinv(reshape(nlp_dtXaAB_dA_sum_xxxx____,[n_x^2,n_x^2]),tolerance_master) * reshape(nlp_dtXaAB_dA_sum_xx__,[n_x^2,1]) , [n_x,n_x] );
clear nlp_dtXaAB_da_sum_xa__ nlp_dtXaAB_da_sum_xaxa____ nlp_dtXaAB_dA_sum_xx__ nlp_dtXaAB_dA_sum_xxxx____ ;
%%%%;
nlp_dtXaAB_sum_pos = 0;
for ni=0:n_i-1;
n_t = n_t_i_(1+ni);
t_t_ = t_it__{1+ni};
X_xt__ = X_ixt___{1+ni};
[ ...
 parameter ...
,nlp_dtXaAB ...
] = ...
PAD_nlp_dtXaAB_0( ...
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
);
nlp_dtXaAB_sum_pos = nlp_dtXaAB_sum_pos + nlp_dtXaAB;
clear n_t t_t_ X_xt_ nlp_dtXaAB ;
end;%for ni=0:n_i-1;
%%%%;
parameter.nlp_dtXaAB_sum_pre = nlp_dtXaAB_sum_pre;
parameter.nlp_dtXaAB_sum_pos = nlp_dtXaAB_sum_pos;
%%%%%%%%;
end;%if ~isempty(strfind(str_update,'a_xa_')) | ~isempty(strfind(str_update,'A_xx_'));

if ~isempty(strfind(str_update,'X_xt__'));
%%%%%%%%;
nlp_tXaABYC_integrated_sum_pre=0;
nlp_tXaABYC_sum_pre = 0;
for ni=0:n_i-1;
n_t = n_t_i_(1+ni);
t_t_ = t_it__{1+ni};
X_xt__ = X_ixt___{1+ni};
ignore_Y_xt__ = ignore_Y_ixt___{1+ni};
Y_xt__ = Y_ixt___{1+ni};
[ ...
 parameter ...
,nlp_tXaABYC_integrated ...
,X_opt_xt__ ...
,nlp_tXaABYC ...
,nlp_dtXaAB ...
,nlp_tXYC ...
,nlp_dtXaAB_dX_xt__ ...
,nlp_dtXaAB_dX_xtxt__ ...
,nlp_tXYC_dX_xt__ ...
,nlp_tXYC_dX_xtxt__ ...
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
nlp_tXaABYC_integrated_sum_pre = nlp_tXaABYC_integrated_sum_pre + nlp_tXaABYC_integrated;
nlp_tXaABYC_sum_pre = nlp_tXaABYC_sum_pre + nlp_tXaABYC;
X_ixt___{1+ni} = X_opt_xt__;
clear n_t t_t_ X_xt__ ignore_Y_xt__ Y_xt__ X_opt_xt__ ;
clear nlp_tXaABYC nlp_tXaABYC_integrated nlp_dtXaAB nlp_tXYC nlp_dtXaAB_dX_xt__ nlp_dtXaAB_dX_xtxt__ nlp_tXYC_dX_xt__ nlp_tXYC_dX_xtxt__ ;
end;%for ni=0:n_i-1;
%%%%;
nlp_tXaABYC_integrated_sum_pos = 0;
nlp_tXaABYC_sum_pos = 0;
for ni=0:n_i-1;
n_t = n_t_i_(1+ni);
t_t_ = t_it__{1+ni};
X_xt__ = X_ixt___{1+ni};
ignore_Y_xt__ = ignore_Y_ixt___{1+ni};
Y_xt__ = Y_ixt___{1+ni};
[ ...
 parameter ...
,nlp_tXaABYC_integrated ...
,X_opt_xt__ ...
,nlp_tXaABYC ...
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
nlp_tXaABYC_integrated_sum_pos = nlp_tXaABYC_integrated_sum_pos + nlp_tXaABYC_integrated;
nlp_tXaABYC_sum_pos = nlp_tXaABYC_sum_pos + nlp_tXaABYC;
clear n_t t_t_ X_xt__ ignore_Y_xt__ Y_xt__ X_opt_xt__ ;
clear nlp_tXaABYC nlp_tXaABYC_integrated ;
end;%for ni=0:n_i-1;
%%%%;
parameter.nlp_tXaABYC_sum_pre = nlp_tXaABYC_sum_pre;
parameter.nlp_tXaABYC_sum_pos = nlp_tXaABYC_sum_pos;
parameter.nlp_tXaABYC_integrated_sum_pre = nlp_tXaABYC_integrated_sum_pre;
parameter.nlp_tXaABYC_integrated_sum_pos = nlp_tXaABYC_integrated_sum_pos;
%%%%%%%%;
end;%if ~isempty(strfind(str_update,'X_xt__'));







