function ...
[ ...
 nlp_tXaABYC_integrated_sum ...
,nlp_tXaABYC_sum ...
,X_opt_ixt___ ...
] = ...
PAD_nlp_itaABYC_strip_0( ...
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

str_thisfunction = 'PAD_nlp_itaABYC_strip_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
nf=0;
rng(0);
n_x = 2;
n_i = 2;
n_t_i_ = zeros(n_i,1);
t_it__ = cell(n_i,1);
ignore_Y_ixt___ = cell(n_i,1);
Y_ixt___ = cell(n_i,1);
for ni=0:n_i-1;
n_t = 8 + ni;
t_t_ = sort(rand(n_t,1));
n_t_i_(1+ni) = n_t;
t_it__{1+ni} = t_t_;
Y_xt__ = randn(n_x,n_t);
Y_ixt___{1+ni} = Y_xt__;
ignore_Y_xt__ = rand(n_x,n_t)< 0.5;
ignore_Y_ixt___{1+ni} = ignore_Y_xt__;
end;%for ni=0:n_i-1;
n_a = 3; a_xa__ = 1*randn(n_x,n_a);
A_xx__ = randn(n_x,n_x);
B_omega = pi/7 ; B_l0 = 1*randn(); B_l1 = 1*randn();
C_omega = pi/5 ; C_l0 = 1*randn(); C_l1 = 1*randn();
X_ixt___ = [];
parameter = struct('type','parameter');
[ ...
 nlp_tXaABYC_integrated_sum ...
,nlp_tXaABYC_sum ...
,X_opt_ixt___ ...
] = ...
PAD_nlp_itaABYC_strip_0( ...
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
disp(sprintf(' %% nlp_tXaABYC_integrated_sum: %0.6f',nlp_tXaABYC_integrated_sum));
disp(sprintf(' %% nlp_tXaABYC_sum: %0.6f',nlp_tXaABYC_sum));
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

if isempty(X_ixt___);
X_ixt___ = cell(n_i,1);
for ni=0:n_i-1;
n_t = n_t_i_(1+ni);
X_ixt___{1+ni} = zeros(n_x,n_t);
end;%for ni=0:n_i-1;
end;%if isempty(X_ixt___);

X_opt_ixt___ = cell(n_i,1);
nlp_tXaABYC_integrated_sum = 0;
nlp_tXaABYC_sum = 0;
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
if isfinite(nlp_tXaABYC_integrated);
nlp_tXaABYC_integrated_sum = nlp_tXaABYC_integrated_sum + nlp_tXaABYC_integrated;
end;%if isfinite(nlp_tXaABYC_integrated);
if isfinite(nlp_tXaABYC);
nlp_tXaABYC_sum = nlp_tXaABYC_sum + nlp_tXaABYC;
end;%if isfinite(nlp_tXaABYC);
X_opt_ixt___{1+ni} = X_opt_xt__;
clear n_t t_t_ X_xt__ ignore_Y_xt__ Y_xt__ X_opt_xt__ ;
clear nlp_tXaABYC nlp_tXaABYC_integrated ;
end;%for ni=0:n_i-1;
