function ...
[ ...
 parameter ...
] = ...
PAD_A_0( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
,Y_xt__ ...
,index_Y_missing_l_ ...
,a0_x_ ...
,a1_x_ ...
,A_xx__ ...
,B_xx__ ...
,C_xx__ ...
);

str_thisfunction = 'PAD_A_0';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_t=[]; end; na=na+1;
if (nargin<1+na); t_t_=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); X_xt__=[]; end; na=na+1;
if (nargin<1+na); Y_xt__=[]; end; na=na+1;
if (nargin<1+na); index_Y_missing_l_=[]; end; na=na+1;
if (nargin<1+na); a0_x_=[]; end; na=na+1;
if (nargin<1+na); a1_x_=[]; end; na=na+1;
if (nargin<1+na); A_xx__=[]; end; na=na+1;
if (nargin<1+na); B_xx__=[]; end; na=na+1;
if (nargin<1+na); C_xx__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
% Calculate negative-log-likelihood for Y (in \Real^{n_x}), ;
% assuming a linear SDE for X: ;
% a(t) = a0_x_ + a1_x_*t ; %<-- linear approximation. ;
% (X(t+dt) - a(t+dt)) - (X(t) - a(t)) = d(X(t)-a(t)) = A * (X(t)-a(t)) * dt + B_inv * dW(t) ; %<-- SDE. ;
% and a simple noise model for Y: ;
% Y(t) = X(t) + C_inv * eps(t) ; %<-- observation noise, with eps ~ N(0,1). ;
%%%%%%%%;

if isempty(t_t_); t_t_ = transpose(0:n_t-1); end;
if isempty(a0_x_); a0_x_ = zeros(n_x,1); end;
if isempty(a1_x_); a1_x_ = zeros(n_x,1); end;
if isempty(X_xt__); X_xt__ = zeros(n_x,n_t); end;
if isempty(Y_xt__); Y_xt__ = X_xt__; end;
if isempty(index_Y_missing_l_); index_Y_missing_l_ = efind(~isfinite(Y_xt__)); end;
if isempty(A_xx__); A_xx__ = zeros(n_x,n_x); end;
if isempty(B_xx__); B_xx__ = ones(n_x,n_x); end;
if isempty(C_xx__); C_xx__ = ones(n_x,n_x); end;

%%%%%%%%;
% given X_xt__ and C_xx__, calculate the negative-log-likelihood for Y_xt__. ;
%%%%%%%%;
det_C = det(C_xx__);
CtCn_xx__ = transpose(C_xx__)*C_xx__;
YX_xt__ = ;

l1 = 3.5; l2 = -1.8;
y = pi/7; c = cos(y); s = sin(y); S__ = [ l1 0 ; 0 l2 ]; U__ = [ +c -s ; +s +c ]; 
B__ = U__ * S__ * transpose(U__) ;
D__ = [ l1*c^2 + l2*s^2 , c*s*(l1-l2) ; c*s*(l1-l2) , l1*s^2 + l2*c^2 ];

n_t = 8;
n_x = 2;
X_xt__ = randn(n_x,n_t);
Y_xt__ = randn(n_x,n_t);
omega = pi/7 ; l0 = randn(); l1 = randn();
[ ...
 ~ ...
,BtBn__ ...
,dw0_BtBn__ ...
,dl0_BtBn__ ...
,dl1_BtBn__ ...
,dw0w0_BtBn__ ...
,dw0l0_BtBn__ ...
,dw0l1_BtBn__ ...
] = ...
PAD_BtBn_0( ...
 [] ...
,omega ...
,l0 ...
,l1 ...
);

nlp = @(X_xt__,Y_xt__,C_xx__) 0.5 * sum( (C_xx__*(Y_xt__ - X_xt__)).^2 , 'all' ) - size(X_xt__,2)*(log(det(C_xx__)) - 0.5*size(X_xt__,1)*log(2*pi));
dwnlp = @(X_xt__,Y_xt__,C_xx__) 0.5 * sum( (C_xx__*(Y_xt__ - X_xt__)).^2 , 'all' ) - size(X_xt__,2)*(log(det(C_xx__)) - 0.5*size(X_xt__,1)*log(2*pi));










if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;







