function ...
[ ...
 parameter ...
,Q_set_i_ ...
,J_set_ ...
,J_off_ ...
,K_set_ ...
,K_off_ ...
] = ...
loop_interact_2( ...
 parameter ...
,A__ ...
,J_set_ ...
,K_set_ ...
);
%%%%%%%%;
% attempts to add and subtract while searching. ;
%%%%%%%%;
str_thisfunction = 'loop_interact_2';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
rng(0);
n_var = 256; n_set = ceil(sqrt(n_var));
A__ = 2*rand(n_var,n_var)-1;
J_set_0_ = transpose(0*n_set + [1:n_set]);
K_set_0_ = transpose(n_var-1*n_set+[1:n_set]);
J_set_1_ = transpose(1*n_set + [1:n_set]);
K_set_1_ = transpose(n_var-2*n_set+[1:n_set]);
A__(J_set_0_,K_set_0_) = -1;
A__(K_set_0_,J_set_0_) = +1;
A__(J_set_1_,K_set_1_) = +1;
A__(K_set_1_,J_set_1_) = -1;
pr_ = randperm(n_var); pc_ = randperm(n_var);
pr_ = 1:n_var; pc_ = 1:n_var;
B__ = A__(pr_,pc_);
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
parameter.flag_check = 1;
parameter.flag_drop = 0;
parameter.flag_lump = 1;
[parameter,Q_set_i_,J_set_,J_off_,K_set_,K_off_] = loop_interact_2(parameter,B__);
n_J_set = numel(J_set_); n_J_off = numel(J_off_);
n_K_set = numel(K_set_); n_K_off = numel(K_off_);
n_gap = 8;
C__ = [ ...
        B__(J_set_,J_set_) , zeros(n_J_set,n_gap) ,   B__(J_set_,K_set_) ; ...
      zeros(n_gap,n_J_set) ,   zeros(n_gap,n_gap) , zeros(n_gap,n_K_set) ; ...
        B__(K_set_,J_set_) , zeros(n_K_set,n_gap) ,   B__(K_set_,K_set_) ; ...
	];
D__ = [ ...
        B__(J_off_,J_off_) , zeros(n_J_off,n_gap) ,   B__(J_off_,K_off_) ; ...
      zeros(n_gap,n_J_off) ,   zeros(n_gap,n_gap) , zeros(n_gap,n_K_off) ; ...
        B__(K_off_,J_off_) , zeros(n_K_off,n_gap) ,   B__(K_off_,K_off_) ; ...
	];
figure(1);clf;figmed;fig80s; Alim_ = [-1,+1];
subplot(1,3,1); imagesc(A__,Alim_); axis image; axisnotick; title('ori');
subplot(1,3,2); imagesc(C__,Alim_); axis image; axisnotick; title('set');
subplot(1,3,3); imagesc(D__,Alim_); axis image; axisnotick; title('off');
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); A__=[]; end; na=na+1;
if (nargin<1+na); J_set_=[]; end; na=na+1;
if (nargin<1+na); K_set_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-6; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'flag_check'); parameter.flag_check = 0; end;
if ~isfield(parameter,'flag_drop'); if  isfield(parameter,'flag_lump'); parameter.flag_drop = 1-parameter.flag_lump; else; parameter.flag_drop = 1; end; end;
if ~isfield(parameter,'flag_lump'); if  isfield(parameter,'flag_drop'); parameter.flag_lump = 1-parameter.flag_drop; else; parameter.flag_lump = 0; end; end;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;
flag_check = parameter.flag_check;
flag_drop = parameter.flag_drop;
flag_lump = parameter.flag_lump;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

[n_row,n_col] = size(A__);
assert(n_row==n_col);

J_all_ = transpose(1:n_row);
if flag_drop; if isempty(J_set_); J_set_ = J_all_; end; end;
if flag_lump; if isempty(J_set_); J_set_ = []; end; end;
J_set_ = J_set_(:);
J_off_ = setdiff(J_all_,J_set_);
K_all_ = transpose(1:n_col);
if flag_drop; if isempty(K_set_); K_set_ = K_all_; end; end;
if flag_lump; if isempty(K_set_); K_set_ = []; end; end;
K_set_ = K_set_(:);
K_off_ = setdiff(K_all_,K_set_);
n_J_set = numel(J_set_); n_J_off = numel(J_off_);
n_K_set = numel(K_set_); n_K_off = numel(K_off_);
flag_continue=1; 
niteration=0;
if  flag_drop; n_iteration = n_J_set + n_K_set; end;%if  flag_drop;
if  flag_lump; n_iteration = n_J_off + n_K_off; end;%if  flag_lump;
Q_set_raw = + sum(A__(J_set_,K_set_),'all') - sum(A__(K_set_,J_set_),'all');
Q_set_upd = Q_set_raw;
J_add_raw_ = + reshape(sum(A__(J_all_,K_set_),2),[n_row,1]) - reshape(sum(A__(K_set_,J_all_),1),[n_row,1]);
J_add_upd_ = J_add_raw_;
K_add_raw_ = + reshape(sum(A__(J_set_,K_all_),1),[n_col,1]) - reshape(sum(A__(K_all_,J_set_),2),[n_col,1]);
K_add_upd_ = K_add_raw_;
Q_set_i_ = zeros(n_iteration,1);
%%%%%%%%%%%%%%%%;
while flag_continue;
%%%%%%%%%%%%%%%%;
if flag_check;
Q_set_raw = + sum(A__(J_set_,K_set_),'all') - sum(A__(K_set_,J_set_),'all');
disp(sprintf(' %% Q_set_raw: %+16.15f Q_set_upd: %+16.15f',Q_set_raw,Q_set_upd));
assert(abs(Q_set_raw-Q_set_upd)<tolerance_master);
end;%if flag_check;
Q_set_nrm = Q_set_upd / max(1,(n_J_set-0)*(n_K_set-0));
Q_set_i_(1+niteration) = Q_set_nrm;
if flag_check;
J_add_raw_ = + reshape(sum(A__(J_all_,K_set_),2),[n_row,1]) - reshape(sum(A__(K_set_,J_all_),1),[n_row,1]);
K_add_raw_ = + reshape(sum(A__(J_set_,K_all_),1),[n_col,1]) - reshape(sum(A__(K_all_,J_set_),2),[n_col,1]);
disp(sprintf(' %% fnorm(J_add_raw_-J_add_upd_): %+16.15f',fnorm(J_add_raw_-J_add_upd_)));
disp(sprintf(' %% fnorm(K_add_raw_-K_add_upd_): %+16.15f',fnorm(K_add_raw_-K_add_upd_)));
assert(fnorm(J_add_raw_-J_add_upd_)<tolerance_master);
assert(fnorm(K_add_raw_-K_add_upd_)<tolerance_master);
end;%if flag_check;
val_J_min_from_set = +Inf; val_K_min_from_set = +Inf;
val_J_max_from_off = -Inf; val_K_max_from_off = -Inf;
[val_J_min_from_set,ij_J_min_from_set] = min((Q_set_upd - J_add_upd_(J_set_)).^2); if isempty(val_J_min_from_set); val_J_min_from_set = +Inf; end;
[val_K_min_from_set,ij_K_min_from_set] = min((Q_set_upd - K_add_upd_(K_set_)).^2); if isempty(val_K_min_from_set); val_K_min_from_set = +Inf; end;
[val_J_max_from_off,ij_J_max_from_off] = max((Q_set_upd + J_add_upd_(J_off_)).^2); if isempty(val_J_max_from_off); val_J_max_from_off = -Inf; end;
[val_K_max_from_off,ij_K_max_from_off] = max((Q_set_upd + K_add_upd_(K_off_)).^2); if isempty(val_K_max_from_off); val_K_max_from_off = -Inf; end;
Q_drop_J_min_from_set = (val_J_min_from_set)/max(1,(n_J_set-1)*(n_K_set-0));
Q_drop_K_min_from_set = (val_K_min_from_set)/max(1,(n_J_set-0)*(n_K_set-1));
Q_fill_J_max_from_off = (val_J_max_from_off)/max(1,(n_J_set+1)*(n_K_set-0));
Q_fill_K_max_from_off = (val_K_max_from_off)/max(1,(n_J_set-0)*(n_K_set+1));
if  flag_drop;
Q_fill_J_max_from_off = +Inf; Q_fill_K_max_from_off = +Inf;
[~,tmp_ij] = min( [Q_drop_J_min_from_set;Q_drop_K_min_from_set;Q_fill_J_max_from_off;Q_fill_K_max_from_off] );
end;%if  flag_drop;
if  flag_lump;
Q_drop_J_min_from_set = -Inf; Q_drop_K_min_from_set = -Inf;
[~,tmp_ij] = max( [Q_drop_J_min_from_set;Q_drop_K_min_from_set;Q_fill_J_max_from_off;Q_fill_K_max_from_off] );
end;%if  flag_lump;
if 0;
elseif tmp_ij==1;
J_upd = J_set_(ij_J_min_from_set);
J_set_ = setdiff(J_set_,J_upd);
J_off_ = [J_off_;J_upd];
Q_set_upd = Q_set_upd - J_add_upd_(J_upd);
K_add_upd_ = K_add_upd_ - transpose(A__(J_upd,:)) + A__(:,J_upd);
n_J_set = n_J_set-1; n_J_off = n_J_off+1;
elseif tmp_ij==2;
K_upd = K_set_(ij_K_min_from_set);
K_set_ = setdiff(K_set_,K_upd);
K_off_ = [K_off_;K_upd];
Q_set_upd = Q_set_upd - K_add_upd_(K_upd);
J_add_upd_ = J_add_upd_ - A__(:,K_upd) + transpose(A__(K_upd,:));
n_K_set = n_K_set-1; n_K_off = n_K_off+1;
elseif tmp_ij==3;
J_upd = J_off_(ij_J_max_from_off);
J_off_ = setdiff(J_off_,J_upd);
J_set_ = [J_set_;J_upd];
Q_set_upd = Q_set_upd + J_add_upd_(J_upd);
K_add_upd_ = K_add_upd_ + transpose(A__(J_upd,:)) - A__(:,J_upd);
n_J_set = n_J_set+1; n_J_off = n_J_off-1;
elseif tmp_ij==4;
K_upd = K_off_(ij_K_max_from_off);
K_off_ = setdiff(K_off_,K_upd);
K_set_ = [K_set_;K_upd];
Q_set_upd = Q_set_upd + K_add_upd_(K_upd);
J_add_upd_ = J_add_upd_ + A__(:,K_upd) - transpose(A__(K_upd,:));
n_K_set = n_K_set+1; n_K_off = n_K_off-1;
end;%if 0;
%%%%;
assert(numel(J_set_)==n_J_set); assert(numel(J_off_)==n_J_off);
assert(numel(K_set_)==n_K_set); assert(numel(K_off_)==n_K_off);
if flag_verbose;
disp(sprintf(' %% niteration %.3d/%.3d: f%d: J_(%.3d,%.3d) K_(%.3d,%.3d) Q [%+0.2f %+0.2f %+0.2f %+0.2f] <-- tmp_ij %d',niteration,n_iteration,flag_drop,n_J_set,n_J_off,n_K_set,n_K_off,[Q_drop_J_min_from_set,Q_drop_K_min_from_set,Q_fill_J_max_from_off,Q_fill_K_max_from_off],tmp_ij));
if (mod(niteration,32)==0);
figure(1);clf;figmed;fig80s; Alim_ = [-1,+1];
subplot(2,2,1); imagesc(A__,Alim_); axis image; axisnotick; title('orig');
subplot(2,2,2); imagesc(A__([J_set_;J_off_],[K_set_;K_off_]),Alim_); axis image; axisnotick; title('reorg');
title(sprintf('niteration %d',niteration));
subplot(2,2,3);
hold on;
plot(J_add_upd_,'kx');
plot(J_set_,J_add_upd_(J_set_),'ro');
hold off;
xlim([1,n_row]); xlabel('row'); title('J_add_upd_','Interpreter','none');
subplot(2,2,4);
hold on;
plot(K_add_upd_,'kx');
plot(K_set_,K_add_upd_(K_set_),'ro');
hold off;
xlim([1,n_col]); xlabel('col'); title('K_add_upd_','Interpreter','none');
sgtitle(sprintf('niteration %d/%d',niteration,n_iteration));
drawnow();
end;%if (mod(niteration,16)==0);
end;%if flag_verbose;
%%%%;
flag_continue = (niteration<n_iteration);
if flag_drop; flag_continue = flag_continue & (n_J_set * n_K_set> 0); end;
if flag_lump; flag_continue = flag_continue & (n_J_off + n_K_off> 0); end;
niteration=niteration+1;
%%%%%%%%%%%%%%%%;
end;%while flag_continue;
%%%%%%%%%%%%%%%%;

if flag_drop; J_off_ = [J_off_;J_set_]; K_off_ = [K_off_;K_set_]; J_set_ = []; K_set_ = []; end;
if flag_lump; J_set_ = [J_set_;J_off_]; K_set_ = [K_set_;K_off_]; J_off_ = []; K_off_ = []; end;

Q_set_i_ = Q_set_i_(1:niteration);

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
