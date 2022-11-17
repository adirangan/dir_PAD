function ...
[ ...
 parameter ...
,t_new_t_ ...
,X_new_xt__ ...
] = ...
PAD_tX_flatten_0( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
);

str_thisfunction = 'PAD_tX_flatten_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
n_t = 8; t_t_ = sort(rand(n_t,1));
t_t_ = [t_t_;t_t_;t_t_]; t_t_ = sort(t_t_); n_t = numel(t_t_);
n_x = 2;
X_xt__ = randn(n_x,n_t);
[ ...
 ~ ...
,t_new_t_ ...
,X_new_xt__ ...
] = ...
PAD_tX_flatten_0( ...
 [] ...
,n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
);
figure(1);clf;figmed;
subplot(1,2,1);plot(t_t_,transpose(X_xt__),'o'); xlabel('time');ylabel('x');
subplot(1,2,2);plot(t_new_t_,transpose(X_new_xt__),'o'); xlabel('time');ylabel('x');
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_t=[]; end; na=na+1;
if (nargin<1+na); t_t_=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); X_xt__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
if ~isfield(parameter,'tolerance_dt'); parameter.tolerance_dt = 1e-12; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
tolerance_master = parameter.tolerance_master;
tolerance_dt = parameter.tolerance_dt;
flag_verbose = parameter.flag_verbose;

if flag_verbose; disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

[~,ij_] = sort(t_t_,'ascend'); index_ = ij_-1;
t_new_t_ = t_t_(1+index_);
X_new_xt__ = X_xt__(:,1+index_);

i_t_ = zeros(n_t,1); ni=0;
for nt=0:n_t-1;
flag_increment = 0;
if (nt==0); end;
if (nt> 0);
tmp_t_pos = t_t_(1+nt+0); tmp_t_pre = t_t_(1+nt-1); tmp_dt = tmp_t_pos - tmp_t_pre;
flag_increment = abs(tmp_dt)>tolerance_dt;
end;%if (nt> 0);
if flag_increment; ni=ni+1; end;
i_t_(1+nt) = ni;
end;%for nt=0:n_t-1;

u_i_ = unique(i_t_); n_u = numel(u_i_);
for nu=0:n_u-1;
tmp_index_ = efind(i_t_==u_i_(1+nu));
n_l = numel(tmp_index_);
if n_l> 0;
tmp_X_new_xt__ = X_new_xt__(:,1+tmp_index_);
X_new_xt__(:,1+tmp_index_) = repmat(mean(tmp_X_new_xt__,2),[1,n_l]);
end;%if n_l> 0;
end;%for nu=0:n_u-1;

if flag_verbose; disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
