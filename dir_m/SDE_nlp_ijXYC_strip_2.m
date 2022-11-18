function ...
[ ...
 nlp_ijXYC ...
,CtCn_xx__ ...
] = ...
SDE_nlp_ijXYC_strip_2( ...
 n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_ixt___ ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_ixj___ ...
,Y_ixj___ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);

if nargin<1;
rng(0);
n_x = 2;
n_i = 3; n_t_i_ = 1024*[1;2;3];
X_ixt___ = cell(n_i,1);
n_j_i_ = zeros(n_i,1);
index_nt_from_nj_i__ = cell(n_i,1);
ignore_Y_ixj___ = cell(n_i,1);
Y_ixj___ = cell(n_i,1);
for ni=0:n_i-1;
n_t = n_t_i_(1+ni);
t_t_ = sort(randn(n_t,1));
X_xt__ = randn(n_x,n_t);
t_it__{1+ni} = t_t_;
X_ixt___{1+ni} = X_xt__;
n_j = ceil(n_t*2);
index_nt_from_nj_ = max(0,min(n_t-1,floor(n_t*rand(n_j,1))));
ignore_Y_xj__ = randn(n_x,n_j)>0.5;
Y_xj__ = randn(n_x,n_j);
n_j_i_(1+ni) = n_j;
index_nt_from_nj_i__{1+ni} = index_nt_from_nj_;
ignore_Y_ixj___{1+ni} = ignore_Y_xj__;
Y_ixj___{1+ni} = Y_xj__;
clear n_t t_t_ X_xt__ n_j index_nt_from_nj_ ignore_Y_xj__ Y_xj__ ;
end;%for ni=0:n_i-1;
C_omega = randn();
C_l0 = randn();
C_l1 = randn();
tmp_t = tic();
[ ...
 nlp_ijXYC ...
,CtCn_xx__ ...
] = ...
SDE_nlp_ijXYC_strip_2( ...
 n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_ixt___ ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_ixj___ ...
,Y_ixj___ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
disp(sprintf(' %% nlp_ijXYC: %0.6f',nlp_ijXYC));
tmp_t = toc(tmp_t); disp(sprintf(' %% SDE_nlp_ijXYC_strip_2: %0.6fs',tmp_t));
disp(sprintf(' %% returning')); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); n_i=[]; end; na=na+1;
if (nargin<1+na); n_t_i_=[]; end; na=na+1;
if (nargin<1+na); t_it__=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); X_ixt___=[]; end; na=na+1;
if (nargin<1+na); n_j_i_=[]; end; na=na+1;
if (nargin<1+na); index_nt_from_nj_i__=[]; end; na=na+1;
if (nargin<1+na); ignore_Y_ixj___=[]; end; na=na+1;
if (nargin<1+na); Y_ixj___=[]; end; na=na+1;
if (nargin<1+na); C_omega=[]; end; na=na+1;
if (nargin<1+na); C_l0=[]; end; na=na+1;
if (nargin<1+na); C_l1=[]; end; na=na+1;

flag_verbose=0;
nlp_ijXYC = 0;
for ni=0:n_i-1;
n_t = n_t_i_(1+ni);
t_t_ = t_it__{1+ni};
X_xt__ = X_ixt___{1+ni};
n_j = n_j_i_(1+ni);
index_nt_from_nj_ = index_nt_from_nj_i__{1+ni};
ignore_Y_xj__ = ignore_Y_ixj___{1+ni};
Y_xj__ = Y_ixj___{1+ni};
[ ...
 nlp_jXYC ...
,CtCn_xx__ ...
] = ...
SDE_nlp_jXYC_strip_2( ...
 n_x ...
,n_t ...
,X_xt__ ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
nlp_ijXYC = nlp_ijXYC + nlp_jXYC;
if flag_verbose; disp(sprintf(' %% ni %d/%d: n_t %d n_j %d nlp_jXYC %+0.6f',ni,n_i,n_t,n_j,nlp_jXYC)); end;
clear n_t t_t_ X_xt__ n_j index_nt_from_nj_ ignore_Y_xj__ Y_xj__ nlp_jXYC;
end;%for ni=0:n_i-1;






