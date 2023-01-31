function ...
[ ...
 parameter ...
] = ...
PAD_fit_0( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_i ...
,n_x ...
,str_var_use_x_ ...
,data_0in_ixt___ ...
,data_0in_measure_ixt___ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
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
if (nargin<1+na); n_t=[]; end; na=na+1;
if (nargin<1+na); t_t_=[]; end; na=na+1;
if (nargin<1+na); n_i=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); str_var_use_x_=[]; end; na=na+1;
if (nargin<1+na); data_0in_ixt___=[]; end; na=na+1;
if (nargin<1+na); data_0in_measure_ixt___=[]; end; na=na+1;
if (nargin<1+na); n_a=[]; end; na=na+1;
if (nargin<1+na); a_xa__=[]; end; na=na+1;
if (nargin<1+na); A_xx__=[]; end; na=na+1;
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
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;
flag_normalize = parameter.flag_normalize;

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
if isempty(n_a); n_a = 3; end;
if isempty(a_xa__); a_xa__=zeros(n_x,n_a); end;
if isempty(A_xx__); A_xx__=zeros(n_x,n_x); end;
if isempty(B_omega); B_omega=0; end;
if isempty(B_l0); B_l0=0; end;
if isempty(B_l1); B_l1=0; end;
if isempty(C_omega); C_omega=0; end;
if isempty(C_l0); C_l0=0; end;
if isempty(C_l1); C_l1=0; end;

data_missing_ixt___ = ~data_0in_measure_ixt___;
index_measure_l_ = efind(data_0in_measure_ixt___);
index_missing_l_ = efind(data_missing_ixt___);
if flag_normalize;
tmp_index_ = efind(data_0in_measure_ixt___);
data_avg = mean(data_0in_ixt___(1+tmp_index_),'all','omitnan');
data_std = std(data_0in_ixt___(1+tmp_index_),1,'all','omitnan');
data_nrm_ixt___ = (data_0in_ixt___ - data_avg)/max(1e-12,data_std);
end;%if flag_normalize;

t_lim_ = [min(t_t_),max(t_t_)];

%%%%%%%%;
% Note that here we only use time-steps where ;
% either one of the variables exists, ;
% or where both exist. ;
% This means that we skip over time-steps with no measured data. ;
% Thus, our first-order approximation will be inaccurate ;
% when the number of skipped time-steps is large. ;
%%%%%%%%;
n_t_use_i_ = zeros(n_i,1);
t_t_use_it__ = cell(n_i,1);
data_use_ixt___ = cell(n_i,1);
data_measure_use_ixt___ = cell(n_i,1);
data_missing_use_ixt___ = cell(n_i,1);
index_measure_il__ = cell(n_i,1);
index_missing_il__ = cell(n_i,1);
for ni=0:n_i-1;
tmp_data_xt__ = squeeze(data_0in_ixt___(1+ni,:,:));
tmp_data_nrm_xt__ = squeeze(data_nrm_ixt___(1+ni,:,:));
tmp_data_measure_xt__ = squeeze(data_0in_measure_ixt___(1+ni,:,:));
index_t_use_ = efind(sum(tmp_data_measure_xt__,1));
n_t_use_i_(1+ni) = numel(index_t_use_);
t_t_use_it__{1+ni} = t_t_(1+index_t_use_);
data_use_ixt___{1+ni} = tmp_data_xt__(:,1+index_t_use_);
data_nrm_use_ixt___{1+ni} = tmp_data_nrm_xt__(:,1+index_t_use_);
data_measure_use_ixt___{1+ni} =  tmp_data_measure_xt__(:,1+index_t_use_);
data_missing_use_ixt___{1+ni} = ~tmp_data_measure_xt__(:,1+index_t_use_);
index_measure_il__{1+ni} = efind(data_measure_use_ixt___{1+ni});
index_missing_il__{1+ni} = efind(data_missing_use_ixt___{1+ni});
end;%for ni=0:n_i-1;

if flag_verbose;
figure(1+nf);nf=nf+1;clf;figsml;
plot(0:n_i-1,n_t_use_i_,'o');
xlabel('ni','Interpreter','none');
ylabel('n_t_use_i_','Interpreter','none');
end;%if flag_verbose;

if flag_verbose;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 8;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
for nx=0:n_x-1;
subplot(1,n_x,1+nx);
hold on;
for ni=0:n_i-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*ni/n_i)));
tmp_n_nage = n_t_use_i_(1+ni);
tmp_t_t_ = t_t_use_it__{1+ni};
tmp_data_xt__ = data_use_ixt___{1+ni};
tmp_data_nrm_xt__ = data_nrm_use_ixt___{1+ni};
tmp_data_measure_xt__ = data_measure_use_ixt___{1+ni};
if tmp_n_nage>2;
tmp_index_ = efind(tmp_data_measure_xt__(1+nx,:));
%plot(tmp_t_t_(1+tmp_index_),tmp_data_xt__(1+nx,1+tmp_index_),'.','MarkerSize',markersize_use,'MarkerFaceColor',c_80s__(1+nc_80s,:));
plot(tmp_t_t_(1+tmp_index_),tmp_data_nrm_xt__(1+nx,1+tmp_index_),'.','MarkerSize',markersize_use,'MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%if tmp_n_nage>2;
end;%for ni=0:n_i-1;
hold off;
xlim(t_lim_); grid on;
xlabel('age','Interpreter','none');
ylabel('value','Interpreter','none');
title(str_var_use_x_{1+nx},'Interpreter','none');
end;%for nx=0:n_x-1;
end;%if flag_verbose;

%%%%%%%%;
% initialization. ;
%%%%%%%%;
ignore_Y_ixt___ = cell(n_i,1);
Y_ixt___ = cell(n_i,1);
X_ixt___ = cell(n_i,1);
t_it__ = cell(n_i,1);
n_t_i_ = zeros(n_i,1);
for ni=0:n_i-1;
Y_ixt___{1+ni} = data_nrm_use_ixt___{1+ni};
ignore_Y_ixt___{1+ni} = zeros(size(Y_ixt___{1+ni}));
X_ixt___{1+ni} = zeros(size(Y_ixt___{1+ni}));
t_it__{1+ni} = t_t_use_it__{1+ni};
n_t_i_(1+ni) = numel(t_it__{1+ni});
end;%for ni=0:n_i-1;
%%%%;
for ni=0:n_i-1;
tmp_index_missing_ = index_missing_il__{1+ni};
Y_ixt___{1+ni}(1+tmp_index_missing_) = X_ixt___{1+ni}(1+tmp_index_missing_);
ignore_Y_ixt___{1+ni}(1+tmp_index_missing_) = 1;
if (flag_verbose>1); disp(sprintf(' %% ni %d/%d: numel missing %d',ni,n_i,numel(tmp_index_missing_))); end;
end;%for ni=0:n_i-1;
%%%%%%%%;

%%%%%%%%;
% fit. ;
%%%%%%%%;
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
PAD_nlp_itXaABYC_update_all_0( ...
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

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

