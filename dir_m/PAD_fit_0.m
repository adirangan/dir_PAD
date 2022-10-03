function ...
[ ...
 parameter ...
] = ...
PAD_fit_0( ...
 parameter ...
,n_age ...
,age_ ...
,n_iid ...
,n_var ...
,data_0in_iva___ ...
,data_0in_measure_iva___ ...
,n_q ...
,q_vq__ ...
,A_vv__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);

str_thisfunction = 'PAD_fit_0';
nf=0;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_age=[]; end; na=na+1;
if (nargin<1+na); age_=[]; end; na=na+1;
if (nargin<1+na); n_iid=[]; end; na=na+1;
if (nargin<1+na); n_var=[]; end; na=na+1;
if (nargin<1+na); data_0in_iva___=[]; end; na=na+1;
if (nargin<1+na); data_0in_measure_iva___=[]; end; na=na+1;
if (nargin<1+na); n_q=[]; end; na=na+1;
if (nargin<1+na); q_vq__=[]; end; na=na+1;
if (nargin<1+na); A_vv__=[]; end; na=na+1;
if (nargin<1+na); B_omega=[]; end; na=na+1;
if (nargin<1+na); B_l0=[]; end; na=na+1;
if (nargin<1+na); B_l1=[]; end; na=na+1;
if (nargin<1+na); C_omega=[]; end; na=na+1;
if (nargin<1+na); C_l0=[]; end; na=na+1;
if (nargin<1+na); C_l1=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'flag_normalize'); parameter.flag_normalize = 1; end;
if ~isfield(parameter,'n_iteration'); parameter.n_iteration = 8; end;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;
flag_normalize = parameter.flag_normalize;
n_iteration = parameter.n_iteration;

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
if isempty(n_q); n_q = 3; end;
if isempty(q_vq__); q_vq__=zeros(n_var,n_q); end;
if isempty(A_vv__); A_vv__=zeros(n_var,n_var); end;
if isempty(B_omega); B_omega=0; end;
if isempty(B_l0); B_l0=0; end;
if isempty(B_l1); B_l1=0; end;
if isempty(C_omega); C_omega=0; end;
if isempty(C_l0); C_l0=0; end;
if isempty(C_l1); C_l1=0; end;

data_missing_iva___ = ~data_0in_measure_iva___;
index_measure_l_ = efind(data_0in_measure_iva___);
index_missing_l_ = efind(data_missing_iva___);
if flag_normalize;
tmp_index_ = efind(data_0in_measure_iva___);
data_avg = mean(data_0in_iva___(1+tmp_index_),'all','omitnan');
data_std = std(data_0in_iva___(1+tmp_index_),1,'all','omitnan');
data_nrm_iva___ = (data_0in_iva___ - data_avg)/max(1e-12,data_std);
end;%if flag_normalize;

age_lim_ = [min(age_),max(age_)];

%%%%%%%%;
% Note that here we only use time-steps where ;
% either one of the variables exists, ;
% or where both exist. ;
% This means that we skip over time-steps with no measured data. ;
% Thus, our first-order approximation will be inaccurate ;
% when the number of skipped time-steps is large. ;
%%%%%%%%;
n_nage_use_i_ = zeros(n_iid,1);
age_use_ia__ = cell(n_iid,1);
data_use_iva___ = cell(n_iid,1);
data_measure_use_iva___ = cell(n_iid,1);
data_missing_use_iva___ = cell(n_iid,1);
index_measure_il__ = cell(n_iid,1);
index_missing_il__ = cell(n_iid,1);
for niid=0:n_iid-1;
tmp_data_va__ = squeeze(data_0in_iva___(1+niid,:,:));
tmp_data_nrm_va__ = squeeze(data_nrm_iva___(1+niid,:,:));
tmp_data_measure_va__ = squeeze(data_0in_measure_iva___(1+niid,:,:));
tmp_nage_use_ = efind(sum(tmp_data_measure_va__,1));
n_nage_use_i_(1+niid) = numel(tmp_nage_use_);
age_use_ia__{1+niid} = age_(1+tmp_nage_use_);
data_use_iva___{1+niid} = tmp_data_va__(:,1+tmp_nage_use_);
data_nrm_use_iva___{1+niid} = tmp_data_nrm_va__(:,1+tmp_nage_use_);
data_use_measure_iva___{1+niid} =  tmp_data_measure_va__(:,1+tmp_nage_use_);
data_use_missing_iva___{1+niid} = ~tmp_data_measure_va__(:,1+tmp_nage_use_);
index_measure_il__{1+niid} = efind(data_use_measure_iva___{1+niid});
index_missing_il__{1+niid} = efind(data_use_missing_iva___{1+niid});
end;%for niid=0:n_iid-1;

if flag_verbose;
figure(1+nf);nf=nf+1;clf;figsml;
plot(0:n_iid-1,n_nage_use_i_,'o');
xlabel('niid','Interpreter','none');
ylabel('n_nage_use','Interpreter','none');
end;%if flag_verbose;

if flag_verbose;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 8;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
for nvar=0:n_var-1;
subplot(1,n_var,1+nvar);
hold on;
for niid=0:n_iid-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*niid/n_iid)));
tmp_n_nage = n_nage_use_i_(1+niid);
tmp_age_ = age_use_ia__{1+niid};
tmp_data_va__ = data_use_iva___{1+niid};
tmp_data_nrm_va__ = data_nrm_use_iva___{1+niid};
tmp_data_measure_va__ = data_use_measure_iva___{1+niid};
if tmp_n_nage>2;
tmp_index_ = efind(tmp_data_measure_va__(1+nvar,:));
%plot(tmp_age_(1+tmp_index_),tmp_data_va__(1+nvar,1+tmp_index_),'.','MarkerSize',markersize_use,'MarkerFaceColor',c_80s__(1+nc_80s,:));
plot(tmp_age_(1+tmp_index_),tmp_data_nrm_va__(1+nvar,1+tmp_index_),'.','MarkerSize',markersize_use,'MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%if tmp_n_nage>2;
end;%for niid=0:n_iid-1;
hold off;
xlim(age_lim_); grid on;
xlabel('age','Interpreter','none');
ylabel('value','Interpreter','none');
end;%for nvar=0:n_var-1;
end;%if flag_verbose;

%%%%%%%%;
% Note that we are using a first-order approximation. ;
% Also, we skip time-steps without meaured data. ;
% If we skip many time-steps in a row, ;
% we will need a more accurate approximation (e.g., second-order). ;
%%%%%%%%;
% Loop: ;
% . Set X_va__ to be zero everywhere. ;
% * Set Y_va__(1+index_missing_) = X_va__(1+index_missing_);
% . Use Y_va__ and X_va__ to estimate C_vv__. ;
% . Use X_va__ and q_vq__ and A_vv__ to estimate B_vv__. ;
% . Use Y_va__ and q_vq__ and A_vv__ and B_vv__ to estimate X_va__. ;
% . Use X_va__ and B_vv__ to estimate q_vq__ and A_vv__. ;
% . Return to step (*). ;
%%%%%%%%;

%%%%%%%%;
% initialization. ;
%%%%%%%%;
Y_iva___ = cell(n_iid,1);
X_iva___ = cell(n_iid,1);
t_ia__ = cell(n_iid,1);
n_t_i_ = zeros(n_iid,1);
for niid=0:n_iid-1;
Y_iva___{1+niid} = data_nrm_use_iva___{1+niid};
X_iva___{1+niid} = zeros(size(Y_iva___{1+niid}));
t_ia__{1+niid} = age_use_ia__{1+niid};
n_t_i_(1+niid) = numel(t_ia__{1+niid});
end;%for niid=0:n_iid-1;

niteration=0; flag_continue = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
while flag_continue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for niid=0:n_iid-1;
tmp_index_missing_ = index_missing_il__{1+niid};
Y_iva___{1+niid}(1+tmp_index_missing_) = X_iva___{1+niid}(1+tmp_index_missing_);
if flag_verbose; disp(sprintf(' %% niid %d/%d: numel missing %d',niid,n_iid,numel(tmp_index_missing_))); end;
end;%for niid=0:n_iid-1;
%%%%%%%%;
n_t_all = sum(n_t_i_);
XY_va__ = zeros(n_var,n_t_all);
na=0;
for niid=0:n_iid-1;
XY_va__(:,1+na+[0:n_t_i_(1+niid)-1]) = Y_iva___{1+niid} - X_iva___{1+niid};
na=na + n_t_i_(1+niid);
end;%for niid=0:n_iid-1;
assert(na==n_t_all);
f_nlp = @(w0l0l1_) PAD_nlp_XYC_strip_0(n_var,n_t_all,XY_va__,w0l0l1_(1+0),w0l0l1_(1+1),w0l0l1_(1+2));
[w0l0l1_opt_,fval_opt] = fminsearch(f_nlp,[C_omega,C_l0,C_l1],optimset('MaxFunEvals',1024));
C_omega = w0l0l1_opt_(1+0); C_l0 = w0l0l1_opt_(1+1); C_l1 = w0l0l1_opt_(1+2);
if flag_verbose;
disp(sprintf(' %% C_omega %+0.2f C_l0 %+0.2f C_l1 %+0.2f',C_omega,C_l0,C_l1));
end;%if flag_verbose;
clear XY_va__;
%%%%%%%%;

n_dt_all = sum(n_t_i_-1);
XY_vda__ = zeros(n_var,n_t_all);
nda=0;
for niid=0:n_iid-1;
[ ...
 nlp ...
,BtBn_vv__ ...
,Q_va__ ...
,P_va__ ...
,dX_vda__ ...
,dQ_vda__ ...
,dP_vda__ ...
,ZX_vda__ ...
,ZQ_vda__ ...
,ZP_vda__ ...
] = ...
PAD_nlp_tXaAB_strip_0( ...
 n_t_i_(1+niid) ...
,t_ia__{1+niid} ...
,n_var ...
,X_iva___{1+niid} ...
,n_q ...
,q_vq__ ...
,A_vv__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
XY_vda__(:,1+nda+[0:n_t_i_(1+niid)-1-1]) = ZP_vda__;
nda=nda + n_t_i_(1+niid)-1;
end;%for niid=0:n_iid-1;
assert(nda==n_dt_all);
f_nlp = @(w0l0l1_) PAD_nlp_XYC_strip_0(n_var,n_dt_all,XY_vda__,w0l0l1_(1+0),w0l0l1_(1+1),w0l0l1_(1+2));
[w0l0l1_opt_,fval_opt] = fminsearch(f_nlp,[B_omega,B_l0,B_l1],optimset('MaxFunEvals',1024));
B_omega = w0l0l1_opt_(1+0); B_l0 = w0l0l1_opt_(1+1); B_l1 = w0l0l1_opt_(1+2);
if flag_verbose;
disp(sprintf(' %% B_omega %+0.2f B_l0 %+0.2f B_l1 %+0.2f',B_omega,B_l0,B_l1));
end;%if flag_verbose;

clear XY_vda__;
%%%%%%%%;
flag_continue = (niteration< n_iteration);
niteration = niteration + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%while flag_continue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

