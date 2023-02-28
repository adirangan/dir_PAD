function ...
[ ...
 nlp_tXaABYC ...
,nlp_dtXaAB ...
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

str_thisfunction = 'PAD_nlp_tXaABYC_strip_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
rng(0);
n_t = 128; t_t_ = sort(rand(n_t,1));
n_x = 2; X_xt__ = randn(n_x,n_t);
Y_xt__ = randn(n_x,n_t); ignore_Y_xt__ = rand(n_x,n_t)< 0.5;
n_a = 2; a_xa__ = randn(n_x,n_a);
A_xx__ = randn(n_x,n_x);
B_omega = pi/7 ; B_l0 = randn(); B_l1 = randn();
C_omega = pi/5 ; C_l0 = randn(); C_l1 = randn();
[ ...
 nlp_tXaABYC ...
,nlp_dtXaAB ...
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
);
disp(sprintf(' %% nlp_tXaABYC: %f',nlp_tXaABYC));
disp(sprintf(' %% nlp_dtXaAB: %f',nlp_dtXaAB));
disp(sprintf(' %% nlp_tXYC: %f',nlp_tXYC));
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

na=0;
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

%%%%%%%%%%%%%%%%;
if (n_t<=1);
nlp_tXaABYC = 0; nlp_dtXaAB = 0; nlp_tXYC = 0;
[ ...
 ~ ...
,BtBn_xx__ ...
] = ...
PAD_BtBn_0( ...
 [] ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
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
end;%if (n_t<=1);
%%%%%%%%%%%%%%%%;
if (n_t> 1);
%%%%%%%%%%%%%%%%;

%%%%;
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
 nlp_tXYC ...
,CtCn_xx__ ...
] = ...
PAD_nlp_tXYC_strip_1( ...
 n_x ...
,n_t ...
,X_xt__ ...
,ignore_Y_xt__ ...
,Y_xt__ ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);

nlp_tXaABYC = nlp_dtXaAB + nlp_tXYC ;

%%%%%%%%%%%%%%%%;
end;%if (n_t> 1);
%%%%%%%%%%%%%%%%;





