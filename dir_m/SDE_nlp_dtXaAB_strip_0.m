function ...
[ ...
 nlp ...
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
SDE_nlp_dtXaAB_strip_0( ...
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

str_thisfunction = 'SDE_nlp_dtXaAB_strip_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
rng(0);
n_t = 8; t_t_ = sort(rand(n_t,1));
n_x = 2; X_xt__ = randn(n_x,n_t);
n_a = 2; a_xa__ = randn(n_x,n_a);
A_xx__ = randn(n_x,n_x);
B_omega = pi/7 ; B_l0 = randn(); B_l1 = randn();
[ ...
 nlp ...
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
SDE_nlp_dtXaAB_strip_0( ...
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
);
disp(sprintf(' %% %f',nlp));
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
if (nargin<1+na); Q_xt__=[]; end; na=na+1;
if (nargin<1+na); P_xt__=[]; end; na=na+1;
if (nargin<1+na); dX_xdt__=[]; end; na=na+1;
if (nargin<1+na); dQ_xdt__=[]; end; na=na+1;
if (nargin<1+na); dP_xdt__=[]; end; na=na+1;
if (nargin<1+na); ZX_xdt__=[]; end; na=na+1;
if (nargin<1+na); ZQ_xdt__=[]; end; na=na+1;
if (nargin<1+na); ZP_xdt__=[]; end; na=na+1;

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
nlp = 0;

%%%%%%%%;
if (n_t<=1);
%<-- do nothing. ;
end;%if (n_t<=1);
%%%%%%%%;
if (n_t> 1);
%%%%;
n_dt = n_t-1; dt_dt_ = diff(t_t_);
if isempty(dX_xdt__); dX_xdt__ = diff(X_xt__,1,2); end;%if isempty(dX_xdt__);
if isempty(Q_xt__);
[ ...
 Q_xt__ ...
] = ...
SDE_Q_ta_strip_0( ...
 n_t ...
,t_t_ ...
,n_x ...
,n_a ...
,a_xa__ ...
);
end;%if isempty(Q_xt__);
if isempty(dQ_xdt__); dQ_xdt__ = diff(Q_xt__,1,2); end;%if isempty(dQ_xdt__);
if isempty(P_xt__); P_xt__ = X_xt__ - Q_xt__; end;%if isempty(P_xt__);
if isempty(dP_xdt__); dP_xdt__ = diff(P_xt__,1,2); end;%if isempty(dP_xdt__);
AX_xdt__ = bsxfun(@times,A_xx__*X_xt__(:,1:end-1),reshape(dt_dt_,[1,n_dt]));
AQ_xdt__ = bsxfun(@times,A_xx__*Q_xt__(:,1:end-1),reshape(dt_dt_,[1,n_dt]));
AP_xdt__ = bsxfun(@times,A_xx__*P_xt__(:,1:end-1),reshape(dt_dt_,[1,n_dt]));
if isempty(ZX_xdt__); ZX_xdt__ = dX_xdt__ - AX_xdt__; end;%if isempty(ZX_xdt__);
if isempty(ZQ_xdt__); ZQ_xdt__ = dQ_xdt__ - AQ_xdt__; end;%if isempty(ZQ_xdt__);
if isempty(ZP_xdt__); ZP_xdt__ = dP_xdt__ - AP_xdt__; end;%if isempty(ZP_xdt__);
tmp_error = fnorm(ZP_xdt__ - (ZX_xdt__ - ZQ_xdt__))/max(1e-12,fnorm(ZP_xdt__));
if (tmp_error>=1e-9); disp(sprintf(' %% Warning, tmp_error %0.16f in %s',tmp_error,str_thisfunction)); end;
%%%%;
%{
nlp = 0;
for ndt=0:n_dt-1;
nlp = nlp + transpose(ZP_xdt__(:,1+ndt)) * ( 0.5 * BtBn_xx__ / dt_dt_(1+ndt) ) * ZP_xdt__(:,1+ndt);
end;%for ndt=0:n_dt-1;
nlp = nlp - n_dt*( (B_l0 + B_l1)/2 - 0.5*n_x*log(2*pi) ) + sum( 0.5*n_x*log(dt_dt_) ) ;
 %}
%nlp = 0.5 * sum( sum( (BtBn_xx__*ZP_xdt__).*ZP_xdt__ , 1 ) ./ reshape(dt_dt_,[1,n_dt]) , 2 ) - n_dt*( (B_l0 + B_l1)/2 - 0.5*n_x*log(2*pi) ) + sum( 0.5*n_x*log(dt_dt_) ) ;
[ ...
 nlp ...
] = ...
SDE_nlp_dtZPB_strip_0( ...
 n_dt ...
,dt_dt_ ...
,n_x ...
,ZP_xdt__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
%%%%;
end;%if (n_t> 1);
%%%%%%%%;






