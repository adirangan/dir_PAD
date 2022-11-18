function ...
[ ...
 parameter ...
,n_t ...
,t_t_ ...
,Q_xt__ ...
,P_xt__ ...
,X_xt__ ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,R_avg ...
,zlim_2x__ ...
] = ...
SDE_generate_data_2( ...
 parameter ...
,n_x ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,BtBn_xx__ ...
,BtBn_inv_xx__ ...
,CtCn_xx__ ...
,CtCn_inv_xx__ ...
,X_ini_x_ ...
);

str_thisfunction = 'SDE_generate_data_2';

if nargin<1;
nf=0;
tolerance_master = 1e-6;
rseed = 1;
n_x = 2; T_ini = 0; T_max = 256; dt_0in = 1/32;
A_xx__ = [-0.3 , +1.0   ; ...
            -1 , -0.5 ] ;
n_a = 3;
a_xa__ = [ +1.00 , -2.00/T_max , +9.50/T_max.^2   ; ...
           -0.50 , +1.25/T_max , -9.25/T_max.^2 ] ;
B_inv__ = [ +0.5 , 0.3  ; ...
	    -0.4 , 0.8] ;
BtBn_inv_xx__ = B_inv__*transpose(B_inv__);
BtBn_xx__ = pinv(BtBn_inv_xx__,tolerance_master);
C_inv__ = [ +0.5 , +0.3  ;   ...
	    -0.5 , +0.7] ;
CtCn_inv_xx__ = C_inv__*transpose(C_inv__);
CtCn_xx__ = pinv(CtCn_inv_xx__,tolerance_master);
%%%%%%%%;
for nflag=0:1;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.dt_avg = dt_0in;
parameter.flag_discrete_vs_exponential = nflag;
parameter.T_ini = T_ini;
parameter.T_max = T_max;
parameter.rseed = rseed;
[ ...
 parameter ...
,n_t ...
,t_t_ ...
,Q_xt__ ...
,P_xt__ ...
,X_xt__ ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,R_avg ...
,zlim_2x__ ...
] = ...
SDE_generate_data_2( ...
 parameter ...
,n_x ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,BtBn_xx__ ...
,BtBn_inv_xx__ ...
,CtCn_xx__ ...
,CtCn_inv_xx__ ...
);
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);plot(sort(diff(t_t_)),'.'); xlabel('nt');ylabel('dt');
subplot(1,2,2);hist(diff(t_t_),128); title('hist');
sgtitle(sprintf('parameter.flag_discrete_vs_exponential = %d',parameter.flag_discrete_vs_exponential),'Interpreter','none');
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
p_row = 2; p_col = 4; np=0;
for pcol=0:p_col-1;
if pcol==0; tmp_t_t_ = t_t_; tmp_W_xt__ = Q_xt__; tmp_str = 'Q_xt__'; end;
if pcol==1; tmp_t_t_ = t_t_; tmp_W_xt__ = P_xt__; tmp_str = 'P_xt__'; end;
if pcol==2; tmp_t_t_ = t_t_; tmp_W_xt__ = X_xt__; tmp_str = 'X_xt__'; end;
if pcol==3; tmp_t_t_ = t_t_(1+index_nt_from_nj_); tmp_W_xt__ = Y_xj__; tmp_str = 'Y_xj__'; end;
subplot(p_row,p_col,1+pcol+0*p_col);
plot(tmp_t_t_,tmp_W_xt__(1,:),'r-',tmp_t_t_,tmp_W_xt__(2,:),'b-');
xlabel('time');ylabel(tmp_str,'Interpreter','none'); xlim([T_ini,T_max]); ylim([min(zlim_2x__(1+0,:)),max(zlim_2x__(1+1,:))])
subplot(p_row,p_col,1+pcol+1*p_col);
s = surfline_0(tmp_W_xt__(1,:),tmp_W_xt__(2,:),tmp_t_t_); set(s,'LineWidth',3);
xlabel(sprintf('%s(1+0,:)',tmp_str),'Interpreter','none'); xlim(zlim_2x__(:,1+0));
ylabel(sprintf('%s(1+1,:)',tmp_str),'Interpreter','none'); ylim(zlim_2x__(:,1+1));
grid on;
end;%for pcol=0:p_col-1;
sgtitle(sprintf('parameter.flag_discrete_vs_exponential = %d',parameter.flag_discrete_vs_exponential),'Interpreter','none');
%%%%%%%%;
end;%for nflag=0:1;
%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); n_a=[]; end; na=na+1;
if (nargin<1+na); a_xa__=[]; end; na=na+1;
if (nargin<1+na); A_xx__=[]; end; na=na+1;
if (nargin<1+na); BtBn_xx__=[]; end; na=na+1;
if (nargin<1+na); BtBn_inv_xx__=[]; end; na=na+1;
if (nargin<1+na); CtCn_xx__=[]; end; na=na+1;
if (nargin<1+na); CtCn_inv_xx__=[]; end; na=na+1;
if (nargin<1+na); X_ini_x_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end;
if (~isfield(parameter,'dt_avg')); parameter.dt_avg = 0.25; end;
if (~isfield(parameter,'flag_discrete_vs_exponential')); parameter.flag_discrete_vs_exponential = 0; end;
if (~isfield(parameter,'T_ini')); parameter.T_ini = 0; end;
if (~isfield(parameter,'T_max')); parameter.T_max = 64; end;
if (~isfield(parameter,'n_j_factor')); parameter.n_j_factor = 2; end;
if (~isfield(parameter,'ignore_factor')); parameter.ignore_factor = 0.5; end;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;
dt_avg = parameter.dt_avg;
flag_discrete_vs_exponential = parameter.flag_discrete_vs_exponential;
T_ini = parameter.T_ini;
T_max = parameter.T_max;
n_j_factor = parameter.n_j_factor;
ignore_factor = parameter.ignore_factor;
rseed = parameter.rseed;
n_t = 1+floor((T_max-T_ini)/dt_avg);
t_t_ = zeros(n_t,1);
rng(rseed);
t_t_(1+0) = T_ini;
for nt=1:n_t-1;
if (flag_discrete_vs_exponential==0); dt=-dt_avg * log(1-rand()); end;
if (flag_discrete_vs_exponential==1); dt=dt_avg; end;
dt_dt_(1+nt-1) = dt;
t_t_(1+nt) = t_t_(1+nt-1) + dt;
end;%for nt=1:n_t-1;
n_dt = n_t-1;
dt_dt_ = diff(t_t_);
n_j = ceil(n_t*n_j_factor);
index_nt_from_nj_ = sort(max(0,min(n_t-1,floor(n_t*rand(n_j,1)))),'ascend');
ignore_Y_xj__ = rand(n_x,n_j)<=ignore_factor;

if flag_verbose; disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

[Psi_xx__,Lambda_xx__] = eig(A_xx__); Psi_inv_xx__ = inv(Psi_xx__);
Lambda_x_ = diag(Lambda_xx__);
%%%%%%%%;
if ( isempty(BtBn_inv_xx__));
BtBn_inv_xx__ = pinv(BtBn_xx__,tolerance_master);
end;%if ( isempty(BtBn_inv_xx__));
B_pinv_xx__ = sqrtm(BtBn_inv_xx__);
if ( isempty(CtCn_inv_xx__));
CtCn_inv_xx__ = pinv(CtCn_xx__,tolerance_master);
end;%if ( isempty(CtCn_inv_xx__));
C_pinv_xx__ = sqrtm(CtCn_inv_xx__);
%%%%%%%%;
if isempty(X_ini_x_); X_ini_x_ = zeros(n_x,1); end;

[ ...
 Q_xt__ ...
] = ...
PAD_Q_ta_strip_0( ...
 n_t ...
,t_t_ ...
,n_x ...
,n_a ...
,a_xa__ ...
);
Q_ini_x_ = Q_xt__(:,1+0);

X_xt__ = zeros(n_x,n_t); X_xt__(:,1+0) = X_ini_x_;
P_xt__ = zeros(n_x,n_t); %<-- P_xt__ = X_xt__ - Q_xt__ ;
P_ini_x_ = X_ini_x_ - Q_ini_x_; P_xt__(:,1+0) = P_ini_x_;
R_dt_ = zeros(n_dt,1);
for nt=1:n_t-1;
t_pre = t_t_(1+nt-1);
t_pos = t_t_(1+nt-0);
ndt = nt-1;
dt = dt_dt_(1+ndt);
X_pre_x_ = X_xt__(:,1+nt-1);
Q_pre_x_ = Q_xt__(:,1+nt-1);
Q_pos_x_ = Q_xt__(:,1+nt-0);
P_pre_x_ = P_xt__(:,1+nt-1);
tmp_exp_Adt_xx__ = Psi_xx__ * diag(exp(+Lambda_x_*dt)) * Psi_inv_xx__;
if flag_discrete_vs_exponential==0;
P_pos_x_ = real( + tmp_exp_Adt_xx__ * ( P_pre_x_ + SDE_sample_int_expnAsBdW_2(dt,Psi_xx__,Lambda_x_,Psi_inv_xx__,B_pinv_xx__,tolerance_master)) );
end;%if flag_discrete_vs_exponential==0;
if flag_discrete_vs_exponential==1;
P_pos_x_ = P_pre_x_ + A_xx__*P_pre_x_*dt + B_pinv_xx__*randn(n_x,1)*sqrt(dt); %<-- discrete ito-step. ;
end;%if flag_discrete_vs_exponential==1;
P_xt__(:,1+nt-0) = P_pos_x_;
X_pos_x_ = P_pos_x_ + Q_pos_x_;
X_xt__(:,1+nt-0) = X_pos_x_;
R_dt_(1+ndt) = fnorm(X_pos_x_ - X_pre_x_)/max(1e-12,fnorm(X_pre_x_));
end;%for nt=1:n_t-1;
R_avg = mean(R_dt_);
%%%%;
[ ...
 n_nj_from_nt_ ...
,index_nj_from_nt__ ...
] = ...
SDE_index_nj_from_nt_0( ...
 n_t ...
,n_j ...
,index_nt_from_nj_ ...
);
Y_xj__ = zeros(n_x,n_j);
for nt=0:n_t-1;
n_nj_from_nt = n_nj_from_nt_(1+nt);
index_nj_from_nt_ = index_nj_from_nt__{1+nt};
if n_nj_from_nt> 0;
Y_xj__(:,1+index_nj_from_nt_) = bsxfun(@plus,X_xt__(:,1+nt),C_pinv_xx__*randn(n_x,n_nj_from_nt));
end;%if n_nj_from_nt> 0;
clear n_nj_from_nt index_nj_from_nt_;
end;%for nt=0:n_t-1;

%%%%%%%%;
% find axis limits. ;
%%%%%%%%;
zlim_2x__ = zeros(2,n_x);
for nx=0:n_x-1;
zlim_2x__(1+0,1+nx) = +Inf;
zlim_2x__(1+0,1+nx) = min(zlim_2x__(1+0,1+nx),min(Q_xt__(1+nx,:),[],'all'));
zlim_2x__(1+0,1+nx) = min(zlim_2x__(1+0,1+nx),min(P_xt__(1+nx,:),[],'all'));
zlim_2x__(1+0,1+nx) = min(zlim_2x__(1+0,1+nx),min(X_xt__(1+nx,:),[],'all'));
zlim_2x__(1+0,1+nx) = min(zlim_2x__(1+0,1+nx),min(Y_xj__(1+nx,:),[],'all'));
zlim_2x__(1+1,1+nx) = -Inf;
zlim_2x__(1+1,1+nx) = max(zlim_2x__(1+1,1+nx),max(Q_xt__(1+nx,:),[],'all'));
zlim_2x__(1+1,1+nx) = max(zlim_2x__(1+1,1+nx),max(P_xt__(1+nx,:),[],'all'));
zlim_2x__(1+1,1+nx) = max(zlim_2x__(1+1,1+nx),max(X_xt__(1+nx,:),[],'all'));
zlim_2x__(1+1,1+nx) = max(zlim_2x__(1+1,1+nx),max(Y_xj__(1+nx,:),[],'all'));
end;%for nx=0:n_x-1;
%%%%;
for nx=0:n_x-1;
zlim_2x__(:,1+nx) = mean(zlim_2x__(:,1+nx)) + 1.25*0.5*diff(zlim_2x__(:,1+nx))*[-1;+1];
end;%for nx=0:n_x-1;
%%%%;

if flag_verbose; disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
