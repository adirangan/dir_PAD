%%%%%%%%;
% setting up test for loop_interact_2. ;
% single antagonistic biclique. ;
%%%%%%%%;
str_thisfunction = 'test_loop_interact_2_0';

str_root = 'data';
if exist('platform.type','file');
tmp_fp = fopen('platform.type'); platform_type = fscanf(tmp_fp,'%s'); fclose(tmp_fp);
if strcmp(platform_type,'access1'); str_root = 'data'; end;
if strcmp(platform_type,'eval1'); str_root = 'home'; end;
end;%if exist('platform.type','file');
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_dolphin',str_root);
dir_mat = sprintf('%s/dir_mat',dir_trunk); if ~exist(dir_mat,'dir'); disp(sprintf(' %% mkdir %d',dir_mat)); mkdir(dir_mat); end;
dir_jpg = sprintf('%s/dir_jpg',dir_trunk); if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %d',dir_jpg)); mkdir(dir_jpg); end;

flag_verbose = 1;
flag_recompute = 0;
flag_replot = 1; nf=0;
tolerance_stderr = 1/64; %<-- colorbar resolution. ;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%n_var_ = ceil(2.^[5:.25:10]); n_n_var = numel(n_var_);
n_var_ = ceil(2.^[12]); n_n_var = numel(n_var_);
M_ = [0.15:0.025:0.65]; n_M = numel(M_);
tmp_str_n = sprintf('n%.2dn%.2d',round(min(log2(n_var_))),round(max(log2(n_var_))));
tmp_str_M = sprintf('M%.2dM%.2d',round(100*min(M_)),round(100*max(M_)));
str_infix = sprintf('%s%s',tmp_str_n,tmp_str_M);
  
%%%%%%%%;
% Here we test flag_drop. ;
%%%%%%%%;
fname_mat = sprintf('%s/test_loop_interact_2_%s.mat',dir_mat,str_infix);
if flag_recompute | ~exist(fname_mat,'file');
disp(sprintf(' %% %s not found, creating',fname_mat));
rseed_ = 0:256; n_rseed = numel(rseed_);
%%%%%%%%;
auc_nMr___ = zeros(n_n_var,n_M,n_rseed);
n_rseed_nM__ = zeros(n_n_var,n_M);
n_l = numel(auc_nMr___); nl=0;
%%%%%%%%;
for nn_var=0:n_n_var-1;
n_var = n_var_(1+nn_var);
for nM=0:n_M-1;
M = M_(1+nM);
tmp_auc_r_ = zeros(n_rseed,1);
tmp_n_rseed = 0;
nrseed=0; flag_continue=1;
while flag_continue;
if (flag_verbose>0); if (mod(nl,1024)==0); disp(sprintf(' %% nl %d/%d',nl,n_l)); end; end;
%%%%;
rseed = rseed_(1+nrseed);
rng(rseed);
n_set_0 = ceil(sqrt((n_var.^2).^M/2));
A__ = 2*(rand(n_var,n_var)>0.5)-1;
J_all_0_ = transpose(1:n_var);
J_set_0_ = transpose(0*n_set_0 + [1:n_set_0]);
J_off_0_ = setdiff(J_all_0_,J_set_0_);
K_all_0_ = transpose(1:n_var);
K_set_0_ = transpose(n_var-1*n_set_0+[1:n_set_0]);
K_off_0_ = setdiff(K_all_0_,K_set_0_);
A__(J_set_0_,K_set_0_) = -1;
A__(K_set_0_,J_set_0_) = +1;
px_n_ = transpose(randperm(n_var)); [~,px_i_] = sort(px_n_);
B__ = A__(px_n_,px_n_); assert(fnorm(A__ - B__(px_i_,px_i_))<1e-6);
%%%%;
parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.flag_check = 0;
parameter.flag_drop = 1;
parameter.flag_lump = 0;
[parameter,Q_set_i_,J_set_,J_off_,K_set_,K_off_] = loop_interact_2(parameter,B__);
J_sol_0_ = px_i_(J_set_0_); K_sol_0_ = px_i_(K_set_0_);
[~,ij_J_sol_from_J_sol_,ij_J_off_from_J_sol_] = intersect(J_sol_0_,J_off_,'stable');
[~,ij_J_sol_from_J_sol_,ij_K_off_from_J_sol_] = intersect(J_sol_0_,K_off_,'stable');
[~,ij_K_sol_from_K_sol_,ij_J_off_from_K_sol_] = intersect(K_sol_0_,J_off_,'stable');
[~,ij_K_sol_from_K_sol_,ij_K_off_from_K_sol_] = intersect(K_sol_0_,K_off_,'stable');
J_cmp_0_ = px_i_(J_off_0_); K_cmp_0_ = px_i_(K_off_0_);
[~,ij_J_cmp_from_J_cmp_,ij_J_off_from_J_cmp_] = intersect(J_cmp_0_,J_off_,'stable');
[~,ij_J_cmp_from_J_cmp_,ij_K_off_from_J_cmp_] = intersect(J_cmp_0_,K_off_,'stable');
[~,ij_K_cmp_from_K_cmp_,ij_J_off_from_K_cmp_] = intersect(K_cmp_0_,J_off_,'stable');
[~,ij_K_cmp_from_K_cmp_,ij_K_off_from_K_cmp_] = intersect(K_cmp_0_,K_off_,'stable');
tmp_auc_JJ = auc_0(ij_J_off_from_J_cmp_,ij_J_off_from_J_sol_);
tmp_auc_KJ = auc_0(ij_K_off_from_J_cmp_,ij_K_off_from_J_sol_);
tmp_auc_JK = auc_0(ij_J_off_from_K_cmp_,ij_J_off_from_K_sol_);
tmp_auc_KK = auc_0(ij_K_off_from_K_cmp_,ij_K_off_from_K_sol_);
tmp_auc = mean([tmp_auc_JJ,tmp_auc_JK,tmp_auc_KJ,tmp_auc_KK]);
%%%%;
tmp_auc_r_(1+nrseed) = tmp_auc;
tmp_n_rseed = tmp_n_rseed + 1;
auc_nMr___(1+nn_var,1+nM,1+nrseed) = tmp_auc;
n_rseed_nM__(1+nn_var,1+nM) = tmp_n_rseed;
tmp_auc_avg = mean(tmp_auc_r_(1:tmp_n_rseed));
tmp_auc_std = std(tmp_auc_r_(1:tmp_n_rseed),1);
flag_stderr_large = (tmp_n_rseed<=4) | (tmp_auc_std/max(1,sqrt(tmp_n_rseed)) >= tolerance_stderr);
flag_continue = (tmp_n_rseed<n_rseed) & (flag_stderr_large);
nrseed = nrseed+1; nl = nl + 1;
clear flag_stderr_large tmp_auc tmp_auc_JJ tmp_auc_JK tmp_auc_KJ tmp_auc_KK;
end;%while flag_continue;
if (flag_verbose>0);
disp(sprintf(' %% nn_var %.2d/%.2d nM %.2d/%.2d <-- tmp_n_rseed %.2d/%.2d <-- tmp_auc_avg %+0.2f tmp_auc_std %+0.2f',nn_var,n_n_var,nM,n_M,tmp_n_rseed,n_rseed,tmp_auc_avg,tmp_auc_std));
end;%if (flag_verbose>0);
clear tmp_auc_r_ tmp_n_rseed tmp_auc_avg tmp_auc_std;
%%%%;
end;%for nM=0:n_M-1;
end;%for nn_var=0:n_n_var-1;
%%%%%%%%;
save(fname_mat ...
,'M_','n_M' ...
,'n_var_','n_n_var' ...
,'rseed_','n_rseed' ...
,'n_l' ...
,'tolerance_stderr' ...
,'auc_nMr___' ...
,'n_rseed_nM__' ...
);
%%%%%%%%;
end;%if ~exist(fname_mat,'file');
if  exist(fname_mat,'file');
tmp_ = load(fname_mat);
end;%if  exist(fname_mat,'file');
%%%%%%%%;
auc_avg_nM__ = sum(tmp_.auc_nMr___,3)./max(1,tmp_.n_rseed_nM__);
auc_std_nM__ = sqrt(max(0,sum(tmp_.auc_nMr___.^2,3)./max(1,tmp_.n_rseed_nM__) - auc_avg_nM__.^2));
auc_p15_nM__ = zeros(tmp_.n_n_var,tmp_.n_M);
auc_p50_nM__ = zeros(tmp_.n_n_var,tmp_.n_M);
auc_p85_nM__ = zeros(tmp_.n_n_var,tmp_.n_M);
for nn_var=0:tmp_.n_n_var-1;
for nM=0:tmp_.n_M-1;
n_rseed = tmp_.n_rseed_nM__(1+nn_var,1+nM);
tmp_auc_r_ = tmp_.auc_nMr___(1+nn_var,1+nM,1:n_rseed);
auc_p15_nM__(1+nn_var,1+nM) = prctile(tmp_auc_r_,15);
auc_p50_nM__(1+nn_var,1+nM) = prctile(tmp_auc_r_,50);
auc_p85_nM__(1+nn_var,1+nM) = prctile(tmp_auc_r_,85);
end;%for nM=0:tmp_.n_M-1;
end;%for nn_var=0:tmp_.n_n_var-1;
%%%%%%%%;

fname_fig_pre = sprintf('%s/test_loop_interact_2_%s_drop_FIGA',dir_jpg,str_infix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
%%%%%%%%;
if n_n_var> 1;
figure(1+nf);nf=nf+1;clf;figsml;colormap('hot');
fontsize_use = 12;
subplot(1,1,1);
imagesc(auc_avg_nM__,[0.5,1.0]); axis image; tmp_c_ = colorbar; set(tmp_c_,'TickLength',[0]);
set(gca,'XTick',1:tmp_.n_M,'XTickLabel',num2str(transpose(tmp_.M_),'%0.3f')); xtickangle(90);
set(gca,'YTick',1:tmp_.n_n_var,'YTickLabel',num2str(transpose(tmp_.n_var_),'%.4d'));
xlabel('M','Interpreter','none');
ylabel('N','Interpreter','none');
set(gca,'TickLength',[0,0]);
set(gca,'FontSize',fontsize_use);
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
%close(gcf);
end;%if n_n_var> 1;
%%%%%%%%;
if n_n_var==1;
figure(1+nf);nf=nf+1;clf;figsml;
c_hot__ = colormap('hot'); n_c_hot = size(c_hot__,1);
markersize_sml = 8; markersize_big = 12; fontsize_use = 12;
subplot(1,1,1);
tmp_auc_avg_M_ = auc_avg_nM__;
tmp_auc_p50_M_ = auc_p50_nM__;
tmp_auc_p15_M_ = auc_p15_nM__;
tmp_auc_p85_M_ = auc_p85_nM__;
hold on;
for nM=0:n_M-1;
tmp_M = tmp_.M_(1+nM);
tmp_auc = tmp_auc_p85_M_(1+nM);
nc_hot = max(0,min(n_c_hot-1,floor(n_c_hot*(tmp_auc-0.5)/(1.0-0.5))));
plot(tmp_M,tmp_auc,'ko','MarkerSize',markersize_sml,'MarkerFaceColor',c_hot__(1+nc_hot,:));
tmp_auc = tmp_auc_p15_M_(1+nM);
nc_hot = max(0,min(n_c_hot-1,floor(n_c_hot*(tmp_auc-0.5)/(1.0-0.5))));
plot(tmp_M,tmp_auc,'ko','MarkerSize',markersize_sml,'MarkerFaceColor',c_hot__(1+nc_hot,:));
tmp_auc = tmp_auc_avg_M_(1+nM);
nc_hot = max(0,min(n_c_hot-1,floor(n_c_hot*(tmp_auc-0.5)/(1.0-0.5))));
plot(tmp_M,tmp_auc,'ko','MarkerSize',markersize_big,'MarkerFaceColor',c_hot__(1+nc_hot,:));
end;%for nM=0:n_M-1;
hold off;
xlim([min(tmp_.M_)-0.025,max(tmp_.M_)+0.025]);
ylim([00.35,1.0]);
xlabel('M'); ylabel('AUC');
set(gca,'XTick',tmp_.M_,'XTickLabel',num2str(transpose(tmp_.M_),'%0.3f')); xtickangle(90);
set(gca,'YTick',00.35:0.05:1.0);
set(gca,'TickLength',[0,0]);
grid on;
title(sprintf('N %d',tmp_.n_var_(1+0)));
set(gca,'FontSize',fontsize_use);
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
%close(gcf);
end;%if n_n_var==1;
%%%%%%%%;
end;%if ~exist(fname_fig_jpg,'file');

for N=[43];%for N=[48,512,1024];
fname_fig_pre = sprintf('%s/test_loop_interact_2_%s_drop_N%.4d_FIGB',dir_jpg,str_infix,N);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figmed;colormap(colormap_nlpvt_twosided);
fontsize_use = 14; Alim = 2.0;
subplot(1,1,1);
rseed=0; n_var = N; M = 0.50;
rng(rseed); n_set_0 = ceil(sqrt((n_var.^2).^M/2));
A__ = 2*(rand(n_var,n_var)>0.5)-1;
J_all_0_ = transpose(1:n_var); J_set_0_ = transpose(0*n_set_0 + [1:n_set_0]); J_off_0_ = setdiff(J_all_0_,J_set_0_);
K_all_0_ = transpose(1:n_var); K_set_0_ = transpose(n_var-1*n_set_0+[1:n_set_0]); K_off_0_ = setdiff(K_all_0_,K_set_0_);
str_x_ = cell(n_var,1); for nvar=0:n_var-1; str_x_{1+nvar} = ' '; end;
str_y_ = cell(n_var,1); for nvar=0:n_var-1; str_y_{1+nvar} = ' '; end;
for nij=transpose(K_set_0_); str_x_{nij} = '-'; end;
for nij=transpose(J_set_0_); str_y_{nij} = '-'; end;
for nij=transpose(J_set_0_); str_x_{nij} = '+'; end;
for nij=transpose(K_set_0_); str_y_{nij} = '+'; end;
A__(J_set_0_,K_set_0_) = -1; A__(K_set_0_,J_set_0_) = +1;
px_n_ = transpose(randperm(n_var)); [~,px_i_] = sort(px_n_);
B__ = A__(px_n_,px_n_); assert(fnorm(A__ - B__(px_i_,px_i_))<1e-6);
subplot(1,2,2);
imagesc(A__,Alim*[-1,+1]);
axis image;
set(gca,'XTick',1:n_var,'XTickLabel',str_x_); xtickangle(90);
set(gca,'YTick',1:n_var,'YTickLabel',str_y_);
set(gca,'TickLength',[0,0]);
title('Rearranged');
set(gca,'FontSize',fontsize_use);
subplot(1,2,1);
imagesc(B__,Alim*[-1,+1]);
axis image;
set(gca,'XTick',1:n_var,'XTickLabel',str_x_(px_n_)); xtickangle(90);
set(gca,'YTick',1:n_var,'YTickLabel',str_y_(px_n_));
set(gca,'TickLength',[0,0]);
title('Original');
set(gca,'FontSize',fontsize_use);
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
%close(gcf);
end;%if ~exist(fname_fig_jpg,'file');
end;%for N=[48,512,1024];

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
if (flag_verbose>0); disp('returning'); end;
return;
