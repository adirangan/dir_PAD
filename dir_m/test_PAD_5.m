%%%%%%%%;
% testing out: ;
% subjects_measurements_quarterly_50th_percentile_Chimpanzee.xlsx ;
%%%%%%%%;

str_species = 'Chimpanzee';

platform_type = 'rusty';
if exist('platform.type','file'); fp=fopen('platform.type'); platform_type = fgetl(fp); fclose(fp); end;
if strcmp(platform_type,'eval1'); string_root = '/home'; end;
if strcmp(platform_type,'access1'); string_root = '/data'; end;
if strcmp(platform_type,'Windows'); string_root = 'C:/Users'; end;

platform_user = 'rangan/dir_bcc';
if exist('platform.user','file'); fp=fopen('platform.user'); platform_user = fgetl(fp); fclose(fp); end;

setup_local;

dir_base = sprintf('%s/%s/dir_PAD',string_root,platform_user);
fname_xls = sprintf('%s/subjects_measurements_quarterly_50th_percentile_%s.xlsx',dir_base,str_species);
flag_replot = 1; nf=0;
flag_recompute = 0;
flag_recollect = 0;

fname_mat = sprintf('%s/subjects_measurements_quarterly_50th_percentile_%s.mat',dir_base,str_species);
%%%%%%%%%%%%%%%%;
if ~exist(fname_mat,'file');
%%%%%%%%%%%%%%%%;

T_ = readtable(fname_xls);

%%%%;
tmp_f = @(c) strcmp(c,'subject');
c_index_subject = min(efind(cellfun(tmp_f,T_.Properties.VariableNames)));
T_subject_ = T_{:,1+c_index_subject};
%%%%;
tmp_f = @(c) strcmp(c,'species');
c_index_species = min(efind(cellfun(tmp_f,T_.Properties.VariableNames)));
T_species_ = T_{:,1+c_index_species};
%%%%;
tmp_f = @(c) strcmp(c,'sex');
c_index_sex = min(efind(cellfun(tmp_f,T_.Properties.VariableNames)));
T_sex_ = T_{:,1+c_index_sex};
%%%%;
tmp_f = @(c) strcmp(c,'social_environment');
c_index_social_environment = min(efind(cellfun(tmp_f,T_.Properties.VariableNames)));
T_social_environment_ = T_{:,1+c_index_social_environment};
%%%%;
tmp_f = @(c) strcmp(c,'housing');
c_index_housing = min(efind(cellfun(tmp_f,T_.Properties.VariableNames)));
T_housing_ = T_{:,1+c_index_housing};
%%%%;
tmp_f = @(c) strcmp(c,'diet');
c_index_diet = min(efind(cellfun(tmp_f,T_.Properties.VariableNames)));
T_diet_ = T_{:,1+c_index_diet};
%%%%;
tmp_f = @(c) strcmp(c,'measurement');
c_index_measurement = min(efind(cellfun(tmp_f,T_.Properties.VariableNames)));
T_measurement_ = T_{:,1+c_index_measurement};
%%%%;

%%%%%%%%;
[ ...
 label_subject_enum_ ...
,n_u_label_subject ...
,u_label_subject_ ...
,index_subject_nu_from_nall_ ...
,n_u_label_subject_ ...
,index_subject_nall_from_nu__ ...
] = ...
label_str_to_enum_1( ...
 T_subject_ ...
);
%%%%%%%%;
[ ...
 label_species_enum_ ...
,n_u_label_species ...
,u_label_species_ ...
,index_species_nu_from_nall_ ...
,n_u_label_species_ ...
,index_species_nall_from_nu__ ...
] = ...
label_str_to_enum_1( ...
 T_species_ ...
);
%%%%%%%%;
[ ...
 label_sex_enum_ ...
,n_u_label_sex ...
,u_label_sex_ ...
,index_sex_nu_from_nall_ ...
,n_u_label_sex_ ...
,index_sex_nall_from_nu__ ...
] = ...
label_str_to_enum_1( ...
 T_sex_ ...
);
%%%%%%%%;
[ ...
 label_social_environment_enum_ ...
,n_u_label_social_environment ...
,u_label_social_environment_ ...
,index_social_environment_nu_from_nall_ ...
,n_u_label_social_environment_ ...
,index_social_environment_nall_from_nu__ ...
] = ...
label_str_to_enum_1( ...
 T_social_environment_ ...
);
%%%%%%%%;
[ ...
 label_housing_enum_ ...
,n_u_label_housing ...
,u_label_housing_ ...
,index_housing_nu_from_nall_ ...
,n_u_label_housing_ ...
,index_housing_nall_from_nu__ ...
] = ...
label_str_to_enum_1( ...
 T_housing_ ...
);
%%%%%%%%;
[ ...
 label_diet_enum_ ...
,n_u_label_diet ...
,u_label_diet_ ...
,index_diet_nu_from_nall_ ...
,n_u_label_diet_ ...
,index_diet_nall_from_nu__ ...
] = ...
label_str_to_enum_1( ...
 T_diet_ ...
);
%%%%%%%%;
[ ...
 label_measurement_enum_ ...
,n_u_label_measurement ...
,u_label_measurement_ ...
,index_measurement_nu_from_nall_ ...
,n_u_label_measurement_ ...
,index_measurement_nall_from_nu__ ...
] = ...
label_str_to_enum_1( ...
 T_measurement_ ...
);

%%%%%%%%;
tmp_f = @(c) strcmp(c,'y0');
c_index_y0 = min(efind(cellfun(tmp_f,T_.Properties.VariableNames)));
c_index_yF = numel(T_.Properties.VariableNames)-1;
n_age = c_index_yF - c_index_y0 + 1;
age_ = 0.25*[0:n_age-1];
%%%%%%%%;
n_iid = n_u_label_subject;
n_var = n_u_label_measurement;
data_iva___ = zeros(n_iid,n_var,n_age);
enum_var_v_ = zeros(n_var,1); str_var_v_ = cell(n_var,1);
enum_iid_i_ = zeros(n_iid,1); str_iid_i_ = cell(n_iid,1);
enum_sex_i_ = zeros(n_iid,1); str_sex_i_ = cell(n_iid,1);
enum_soc_i_ = zeros(n_iid,1); str_soc_i_ = cell(n_iid,1);
enum_hou_i_ = zeros(n_iid,1); str_hou_i_ = cell(n_iid,1);
enum_die_i_ = zeros(n_iid,1); str_die_i_ = cell(n_iid,1);
%%%%%%%%;
for nvar=0:n_var-1;
index_var_ = index_measurement_nall_from_nu__{1+nvar};
enum_var_v_(1+nvar) = label_measurement_enum_(1+index_var_(1+0)); 
str_var_v_{1+nvar} = T_measurement_{1+index_var_(1+0)};
end;%for nvar=0:n_var-1;
%%%%%%%%;
for niid=0:n_iid-1;
index_iid_ = index_subject_nall_from_nu__{1+niid};
enum_iid_i_(1+niid) = label_subject_enum_(1+index_iid_(1+0)); 
str_iid_i_{1+niid} = T_subject_{1+index_iid_(1+0)};
enum_sex_i_(1+niid) = label_sex_enum_(1+index_iid_(1+0)); 
str_sex_i_{1+niid} = T_sex_{1+index_iid_(1+0)};
enum_soc_i_(1+niid) = label_social_environment_enum_(1+index_iid_(1+0)); 
str_soc_i_{1+niid} = T_social_environment_{1+index_iid_(1+0)};
enum_hou_i_(1+niid) = label_housing_enum_(1+index_iid_(1+0)); 
str_hou_i_{1+niid} = T_housing_{1+index_iid_(1+0)};
enum_die_i_(1+niid) = label_diet_enum_(1+index_iid_(1+0)); 
str_die_i_{1+niid} = T_diet_{1+index_iid_(1+0)};
end;%for niid=0:n_iid-1;
%%%%%%%%;
flag_verbose=1;
for niid=0:n_iid-1;
index_iid_ = index_subject_nall_from_nu__{1+niid};
if (flag_verbose); if (mod(niid,16)==0); disp(sprintf(' %% niid %d/%d',niid,n_iid)); end; end;
for nvar=0:n_var-1;
index_var_ = index_measurement_nall_from_nu__{1+nvar};
tmp_index_ = intersect(index_iid_,index_var_);
if (numel(tmp_index_)>1);
disp(sprintf(' %% Warning, multiple etries for niid %d nvar %d',niid,nvar));
end;%if (numel(tmp_index_)>1);
if (numel(tmp_index_)>=1);
tmp_index = tmp_index_(1+0);
for nage=0:n_age-1;
ny = c_index_y0 + nage;
data_iva___(1+niid,1+nvar,1+nage) = T_{1+tmp_index,1+ny};
end;%for nage=0:n_age-1;
end;%if (numel(tmp_index_)>=1);
end;%for nvar=0:n_var-1;
end;%for niid=0:n_iid-1;
%%%%;
year_ = age_; n_y = n_age; n_year = n_age;
%%%%%%%%;
save( fname_mat ...
      ,'T_diet_' ...
      ,'T_housing_' ...
      ,'T_measurement_' ...
      ,'T_sex_' ...
      ,'T_social_environment_' ...
      ,'T_species_' ...
      ,'T_subject_' ...
      ,'age_' ...
      ,'c_index_diet' ...
      ,'c_index_housing' ...
      ,'c_index_measurement' ...
      ,'c_index_sex' ...
      ,'c_index_social_environment' ...
      ,'c_index_species' ...
      ,'c_index_subject' ...
      ,'c_index_y0' ...
      ,'c_index_yF' ...
      ,'data_iva___' ...
      ,'dir_base' ...
      ,'enum_die_i_' ...
      ,'enum_hou_i_' ...
      ,'enum_iid_i_' ...
      ,'enum_sex_i_' ...
      ,'enum_soc_i_' ...
      ,'enum_var_v_' ...
      ,'flag_verbose' ...
      ,'fname_xls' ...
      ,'index_diet_nall_from_nu__' ...
      ,'index_diet_nu_from_nall_' ...
      ,'index_housing_nall_from_nu__' ...
      ,'index_housing_nu_from_nall_' ...
      ,'index_iid_' ...
      ,'index_measurement_nall_from_nu__' ...
      ,'index_measurement_nu_from_nall_' ...
      ,'index_sex_nall_from_nu__' ...
      ,'index_sex_nu_from_nall_' ...
      ,'index_social_environment_nall_from_nu__' ...
      ,'index_social_environment_nu_from_nall_' ...
      ,'index_species_nall_from_nu__' ...
      ,'index_species_nu_from_nall_' ...
      ,'index_subject_nall_from_nu__' ...
      ,'index_subject_nu_from_nall_' ...
      ,'index_var_' ...
      ,'label_diet_enum_' ...
      ,'label_housing_enum_' ...
      ,'label_measurement_enum_' ...
      ,'label_sex_enum_' ...
      ,'label_social_environment_enum_' ...
      ,'label_species_enum_' ...
      ,'label_subject_enum_' ...
      ,'n_age' ...
      ,'n_iid' ...
      ,'n_u_label_diet' ...
      ,'n_u_label_diet_' ...
      ,'n_u_label_housing' ...
      ,'n_u_label_housing_' ...
      ,'n_u_label_measurement' ...
      ,'n_u_label_measurement_' ...
      ,'n_u_label_sex' ...
      ,'n_u_label_sex_' ...
      ,'n_u_label_social_environment' ...
      ,'n_u_label_social_environment_' ...
      ,'n_u_label_species' ...
      ,'n_u_label_species_' ...
      ,'n_u_label_subject' ...
      ,'n_u_label_subject_' ...
      ,'n_var' ...
      ,'n_y' ...
      ,'n_year' ...
      ,'year_' ...
      ,'str_die_i_' ...
      ,'str_hou_i_' ...
      ,'str_iid_i_' ...
      ,'str_sex_i_' ...
      ,'str_soc_i_' ...
      ,'str_var_v_' ...
      ,'u_label_diet_' ...
      ,'u_label_housing_' ...
      ,'u_label_measurement_' ...
      ,'u_label_sex_' ...
      ,'u_label_social_environment_' ...
      ,'u_label_species_' ...
      ,'u_label_subject_' ...
      );

%%%%%%%%%%%%%%%%;
end;%if ~exist(fname_mat,'file');
%%%%%%%%%%%%%%%%;
if  exist(fname_mat,'file');
load(fname_mat);
end;%if  exist(fname_mat,'file');

%%%%%%%%;
% reload again in case string_root was overwritten by fname_mat. ;
%%%%%%%%;
platform_type = 'rusty';
if exist('platform.type','file'); fp=fopen('platform.type'); platform_type = fgetl(fp); fclose(fp); end;
if strcmp(platform_type,'eval1'); string_root = '/home'; end;
if strcmp(platform_type,'access1'); string_root = '/data'; end;
if strcmp(platform_type,'Windows'); string_root = 'C:/Users'; end;
platform_user = 'rangan/dir_bcc';
if exist('platform.user','file'); fp=fopen('platform.user'); platform_user = fgetl(fp); fclose(fp); end;
setup_local;
dir_base = sprintf('%s/%s/dir_PAD',string_root,platform_user);
  
data_measure_iva___ = isfinite(data_iva___) & (data_iva___~=0);
data_measure_i_ = squeeze(sum(sum(data_measure_iva___,3),2));
data_measure_v_ = squeeze(sum(sum(data_measure_iva___,3),1));
data_measure_a_ = squeeze(sum(sum(data_measure_iva___,2),1));
%%%%%%%%;
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(2,2,1);
plot(age_,data_measure_a_,'o');
xlabel('age');ylabel('data_measure_a_','Interpreter','none');
subplot(2,2,[3,4]);
plot(1:n_var,data_measure_v_,'o');
xlabel('var');ylabel('data_measure_v_','Interpreter','none');
set(gca,'XTick',1:n_var,'XTickLabel',u_label_measurement_); xtickangle(90);
subplot(2,2,2);
plot(1:n_iid,data_measure_i_,'o');
xlabel('iid');ylabel('data_measure_i_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

age_lim_ = [min(age_),max(age_)];
n_var_use = 2;
%str_var_use_v_ = {'Body Weight','Glucose'};
str_var_use_v_ = {'Red Blood Cells','Mean Corpuscular Hemoglobin'};
nv0 = efind(strcmp(u_label_measurement_,str_var_use_v_{1+0}));
nv1 = efind(strcmp(u_label_measurement_,str_var_use_v_{1+1}));
%%%%%%%%;
% Note that nv0 increases roughly linearly until around age 10-15 or so. ;
% You can try: ;
% plot(age_,squeeze(data_iva___(:,1+[nv0],:)),'o'); grid on;
% This means that a longevity study might start afterwards (e.g., age 15+). ;
%%%%%%%%;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Now plot the data for a single individual (niid==143). ;
%%%%%%%%;
niid = 143;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
plot(age_,squeeze(data_iva___(1+niid,1+[nv0],:)),'ko');
plot(age_,squeeze(data_iva___(1+niid,1+[nv1],:)),'ro');
hold off;
xlabel('age'); ylabel('value');
xlim(age_lim_); grid on;
legend({u_label_measurement_{1+nv0},u_label_measurement_{1+nv1}},'Location','SouthEast');
end;%if flag_disp;

n_shuffle = 64; %<-- increase later. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now step through all variable-pairs. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nv0=0:n_var-1; for nv1=nv0+1:n_var-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose); if ( ((nv0==n_var-2) & (nv1==n_var-1)) | mod(nv1+nv0*n_var,16)==0 ); disp(sprintf(' %% nv0 nv1 (%.2d,%.2d)',nv0,nv1)); end; end;

%%%%%%%%;
% Now extract times and measurements for all the individuals. ;
%%%%%%%%;
n_x = 2;
n_i = n_iid;
n_t_i_ = zeros(n_i,1);
t_it__ = cell(n_i,1);
n_j_i_ = zeros(n_i,1);
index_nt_from_nj_i__ = cell(n_i,1);
ignore_Y_tru_ixj___ = cell(n_i,1);
Y_tru_ixj___ = cell(n_i,1);
for ni=0:n_i-1;
tmp_retain_Y_xt__ = squeeze(data_measure_iva___(1+ni,1+[nv0,nv1],:));
tmp_ignore_Y_xt__ = ~tmp_retain_Y_xt__;
tmp_Y_xt__= squeeze(data_iva___(1+ni,1+[nv0,nv1],:));
tmp_retain_Y_t_ = reshape(sum(tmp_retain_Y_xt__,1),[n_age,1]);
tmp_index_ = efind(tmp_retain_Y_t_> 0);
t_t_ = age_(1+tmp_index_); n_t = numel(t_t_);
n_j = n_t; %<-- no repeated measurements for PAD data. ;
index_nt_from_nj_ = transpose(0:n_t-1); %<-- no repeated measurements for PAD data. ;
ignore_Y_tru_xj__ = tmp_ignore_Y_xt__(:,1+tmp_index_);
Y_tru_xj__ = tmp_Y_xt__(:,1+tmp_index_);
index_nt_from_nj_i__{1+ni} = index_nt_from_nj_;
ignore_Y_tru_ixj___{1+ni} = ignore_Y_tru_xj__;
Y_tru_ixj___{1+ni} = Y_tru_xj__;
n_j_i_(1+ni) = n_j;
n_t_i_(1+ni) = n_t;
t_it__{1+ni} = t_t_;
clear tmp_retain_Y_xt__ tmp_ignore_Y_xt__ tmp_Y_xt__ tmp_retain_Y_t_ tmp_index_ t_t_ n_j index_nt_from_nj_ ignore_Y_tru_xj__ Y_tru_xj__ ;
end;%for ni=0:n_i-1;
%%%%%%%%;
% remove monkeys with fewer than 3 measurements. ;
%%%%%%%%;
tmp_index_ = efind(n_j_i_>=3);
n_i = numel(tmp_index_);
n_t_i_ = n_t_i_(1+tmp_index_);
t_it__ = t_it__(1+tmp_index_);
n_j_i_ = n_j_i_(1+tmp_index_);
index_nt_from_nj_i__ = index_nt_from_nj_i__(1+tmp_index_);
ignore_Y_tru_ixj___ = ignore_Y_tru_ixj___(1+tmp_index_);
Y_tru_ixj___ = Y_tru_ixj___(1+tmp_index_);
%%%%%%%%;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Now plot just to check. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
n_iteration_impute = 8;
%%%%%%%%;
for ni=0:n_i-1;
t_t_ = t_it__{1+ni};
index_nt_from_nj_ = index_nt_from_nj_i__{1+ni};
Y_tru_xj__ = Y_tru_ixj___{1+ni};
ignore_Y_tru_xj__ = ignore_Y_tru_ixj___{1+ni};
retain_Y_tru_xj__ = ~ignore_Y_tru_xj__;
t_j_ = t_t_(1+index_nt_from_nj_);
tmp_index_retain_ = efind(retain_Y_tru_xj__);
tmp_index_retain_0_ = efind(retain_Y_tru_xj__(1+0,:));
tmp_index_retain_1_ = efind(retain_Y_tru_xj__(1+1,:));
tmp_index_ignore_ = efind(ignore_Y_tru_xj__);
tmp_index_ignore_0_ = efind(ignore_Y_tru_xj__(1+0,:));
tmp_index_ignore_1_ = efind(ignore_Y_tru_xj__(1+1,:));
Y_tru_xj__(1+tmp_index_ignore_) = 0;
for niteration_impute=0:n_iteration_impute-1;
[tmp_U__,tmp_S__,tmp_V__] = svds(Y_tru_xj__,1);
tmp_USV__ = tmp_U__*tmp_S__*transpose(tmp_V__);
Y_tru_xj__(1+tmp_index_ignore_) = tmp_USV__(1+tmp_index_ignore_);
end;%for niteration_impute=0:n_iteration_impute-1;
%%%%;
subplot(1,3,1);
hold on;
plot(t_j_(1+tmp_index_retain_0_),Y_tru_xj__(1+0,1+tmp_index_retain_0_),'k.');
plot(t_j_(1+tmp_index_ignore_0_),Y_tru_xj__(1+0,1+tmp_index_ignore_0_),'c.');
hold off;
xlabel('age'); xlim(age_lim_);
ylabel(sprintf('%s',u_label_measurement_{1+nv0}),'Interpreter','none');
%%%%;
subplot(1,3,2);
hold on;
plot(t_j_(1+tmp_index_retain_1_),Y_tru_xj__(1+1,1+tmp_index_retain_1_),'k.');
plot(t_j_(1+tmp_index_ignore_1_),Y_tru_xj__(1+1,1+tmp_index_ignore_1_),'g.');
hold off;
xlabel('age'); xlim(age_lim_);
ylabel(sprintf('%s',u_label_measurement_{1+nv1}),'Interpreter','none');
%%%%;
subplot(1,3,3);
hold on;
tmp_index_ = intersect(tmp_index_retain_0_,tmp_index_retain_1_);
plot(Y_tru_xj__(1+0,1+tmp_index_),Y_tru_xj__(1+1,1+tmp_index_),'k.');
tmp_index_ = intersect(tmp_index_retain_0_,tmp_index_ignore_1_);
plot(Y_tru_xj__(1+0,1+tmp_index_),Y_tru_xj__(1+1,1+tmp_index_),'g.');
tmp_index_ = intersect(tmp_index_ignore_0_,tmp_index_retain_1_);
plot(Y_tru_xj__(1+0,1+tmp_index_),Y_tru_xj__(1+1,1+tmp_index_),'c.');
tmp_index_ = intersect(tmp_index_ignore_0_,tmp_index_ignore_1_);
plot(Y_tru_xj__(1+0,1+tmp_index_),Y_tru_xj__(1+1,1+tmp_index_),'r.');
hold off;
xlabel(sprintf('%s',u_label_measurement_{1+nv0}),'Interpreter','none');
ylabel(sprintf('%s',u_label_measurement_{1+nv1}),'Interpreter','none');
%%%%;
if (ni==n_i-1 | (mod(ni,32)==0)); sgtitle(sprintf('ni %d/%d',ni,n_i),'Interpreter','none'); drawnow(); end;
%%%%%%%%;
end;%for ni=0:n_i-1;
%%%%%%%%;
end;%if flag_disp;

%%%%%%%%;
% Make some extra directories to store the results. ;
%%%%%%%%;
dir_PAD_infix_0 = sprintf('%s',str_species);
dir_PAD_infix_1 = sprintf('%s_nv%.2d',str_species,nv0);
dir_PAD_infix_2 = sprintf('%s_nv%.2dnv%.2d',str_species,nv0,nv1);
%%%%;
dir_PAD_output_0_mat = sprintf('%s/dir_PAD_%s_mat',dir_base,dir_PAD_infix_0);
dir_PAD_output_0_jpg = sprintf('%s/dir_PAD_%s_jpg',dir_base,dir_PAD_infix_0);
if ~exist(dir_PAD_output_0_mat,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_0_mat)); mkdir(dir_PAD_output_0_mat); end;
if ~exist(dir_PAD_output_0_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_0_jpg)); mkdir(dir_PAD_output_0_jpg); end;
%%%%;
dir_PAD_output_1_mat = sprintf('%s/dir_PAD_%s_mat',dir_PAD_output_0_mat,dir_PAD_infix_1);
dir_PAD_output_1_jpg = sprintf('%s/dir_PAD_%s_jpg',dir_PAD_output_0_jpg,dir_PAD_infix_1);
if ~exist(dir_PAD_output_1_mat,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_1_mat)); mkdir(dir_PAD_output_1_mat); end;
if ~exist(dir_PAD_output_1_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_1_jpg)); mkdir(dir_PAD_output_1_jpg); end;
%%%%;
dir_PAD_output_2_mat = sprintf('%s/dir_PAD_%s_mat',dir_PAD_output_1_mat,dir_PAD_infix_2);
dir_PAD_output_2_jpg = sprintf('%s/dir_PAD_%s_jpg',dir_PAD_output_1_jpg,dir_PAD_infix_2);
if ~exist(dir_PAD_output_2_mat,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_2_mat)); mkdir(dir_PAD_output_2_mat); end;
if ~exist(dir_PAD_output_2_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_2_jpg)); mkdir(dir_PAD_output_2_jpg); end;
%%%%;

for nshuffle=0:n_shuffle-1+1;
str_nshuffle = sprintf('s%.4d',nshuffle);
tmp_fname_pre = sprintf('%s/%s_%s',dir_PAD_output_2_mat,dir_PAD_infix_2,str_nshuffle);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre);
%%%%%%%%;
if ~tmp_flag_skip;
%%%%%%%%;
tmp_t = tic();
tolerance_master = 1e-2;
n_a = 2;
flag_regularize_eccentricity_simultaneous = 1;
n_iteration_BtBn = 32;
MaxFunEvals_use_BtBn = 128;
MaxFunEvals_use_simultaneous = 0;
%%%%%%%%;
% recover SDE parameters from data. ;
%%%%%%%%;
X_est_ixt___=[];
a_est_xa__=zeros(n_x,n_a);
A_est_xx__=zeros(n_x,n_x);
B_est_omega=0;
B_est_l0=-10;
B_est_l1=-10;
C_est_omega=0;
C_est_l0=-10;
C_est_l1=-10;
%%%%;
parameter_est = struct('type','parameter_est');
parameter_est.tolerance_master = tolerance_master;
parameter_est.flag_verbose = flag_verbose-0;
parameter_est.flag_disp = flag_verbose-1;
parameter_est.flag_normalize_Y = 1;
parameter_est.nshuffle_Y = nshuffle;
%parameter_est.flag_regularize_eccentricity_BtBn = 1;
parameter_est.flag_regularize_eccentricity_simultaneous = flag_regularize_eccentricity_simultaneous;
parameter_est.n_iteration_BtBn = n_iteration_BtBn;
parameter_est.MaxFunEvals_use_BTBn = MaxFunEvals_use_BtBn;
parameter_est.MaxFunEvals_use_simultaneous = MaxFunEvals_use_simultaneous;
[ ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_est_ixt___ ...
,n_a ...
,a_est_xa__ ...
,A_est_xx__ ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_nrm_prm_ixj___ ...
,Y_tru_nrm_prm_ixj___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
] = ...
SDE_nlp_ijXaABYC_update_all_1( ...
 parameter_est ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_est_ixt___ ...
,n_a ...
,a_est_xa__ ...
,A_est_xx__ ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);
%%%%%%%%;
[ ...
 ~ ...
,BtBn_est_xx__ ...
] = ...
SDE_BtBn_0( ...
 [] ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
);
BtBn_est_inv_xx__ = pinv(BtBn_est_xx__,tolerance_master);
%%%%;
[ ...
 ~ ...
,CtCn_est_xx__ ...
] = ...
SDE_BtBn_0( ...
 [] ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);
CtCn_est_inv_xx__ = pinv(CtCn_est_xx__,tolerance_master);
%%%%;
save(tmp_fname_mat ...
     ,'nv0','nv1' ...
     ,'nshuffle','str_nshuffle','tmp_fname_pre' ...
     ,'str_species' ...
     ,'dir_PAD_infix_0','dir_PAD_output_0_mat','dir_PAD_output_0_jpg' ...
     ,'dir_PAD_infix_1','dir_PAD_output_1_mat','dir_PAD_output_1_jpg' ...
     ,'dir_PAD_infix_2','dir_PAD_output_2_mat','dir_PAD_output_2_jpg' ...
     ,'tolerance_master','n_a','flag_regularize_eccentricity_simultaneous','n_iteration_BtBn','MaxFunEvals_use_BtBn','MaxFunEvals_use_simultaneous' ...
     ,'parameter_est','n_x','n_a','a_est_xa__','A_est_xx__','B_est_omega','B_est_l0','B_est_l1','n_j_i_','C_est_omega','C_est_l0','C_est_l1' ...
     ,'BtBn_est_xx__','CtCn_est_xx__' ...
     );
close_fname_tmp(tmp_fname_pre);
%%%%%%%%;
if nshuffle==0;
%%%%%%%%;
% Now visualize the results. ;
%%%%%%%%;
tmp_fname_fig_pre = sprintf('%s/%s_%s_FIGA',dir_PAD_output_2_jpg,dir_PAD_infix_2,str_nshuffle);
tmp_fname_fig_jpg = sprintf('%s.jpg',tmp_fname_fig_pre);
if flag_replot | ~exist(tmp_fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;
%%%%;
subplot(1,2,1);
hold on;
for ni=0:n_i-1;
t_t_ = t_it__{1+ni};
X_est_xt__ = X_est_ixt___{1+ni};
Y_tru_nrm_prm_xj__ = Y_tru_nrm_prm_ixj___{1+ni};
ignore_Y_tru_nrm_prm_xj__ = ignore_Y_tru_nrm_prm_ixj___{1+ni};
tmp_index_ = efind(~ignore_Y_tru_nrm_prm_xj__(1+0,:));
plot(t_t_(1+tmp_index_),Y_tru_nrm_prm_xj__(1+0,1+tmp_index_),'ko');
plot(t_t_,X_est_xt__(1+0,:),'rx');
end;%for ni=0:n_i-1;
hold off;
xlim(age_lim_);xlabel('time');
ylabel(u_label_measurement_{1+nv0},'Interpreter','none');
%%%%;
subplot(1,2,2);
hold on;
for ni=0:n_i-1;
t_t_ = t_it__{1+ni};
X_est_xt__ = X_est_ixt___{1+ni};
Y_tru_nrm_prm_xj__ = Y_tru_nrm_prm_ixj___{1+ni};
ignore_Y_tru_nrm_prm_xj__ = ignore_Y_tru_nrm_prm_ixj___{1+ni};
tmp_index_ = efind(~ignore_Y_tru_nrm_prm_xj__(1+1,:));
plot(t_t_(1+tmp_index_),Y_tru_nrm_prm_xj__(1+1,1+tmp_index_),'ko');
plot(t_t_,X_est_xt__(1+1,:),'rx');
end;%for ni=0:n_i-1;
hold off;
xlim(age_lim_);xlabel('time');
ylabel(u_label_measurement_{1+nv1},'Interpreter','none');
%%%%;
sgtitle(tmp_fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',tmp_fname_fig_pre));
print('-djpeg',tmp_fname_fig_jpg);
close(gcf);
%%%%;
end;%if ~exist(tmp_fname_fig_jpg,'file');
%%%%%%%%;
tmp_fname_fig_pre = sprintf('%s/%s_%s_FIGB',dir_PAD_output_2_jpg,dir_PAD_infix_2,str_nshuffle);
tmp_fname_fig_jpg = sprintf('%s.jpg',tmp_fname_fig_pre);
if flag_replot | ~exist(tmp_fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;
%%%%;
subplot(1,1,1);
hold on;
for ni=0:n_i-1;
t_t_ = t_it__{1+ni};
X_est_xt__ = X_est_ixt___{1+ni};
Y_tru_nrm_prm_xj__ = Y_tru_nrm_prm_ixj___{1+ni};
ignore_Y_tru_nrm_prm_xj__ = ignore_Y_tru_nrm_prm_ixj___{1+ni};
tmp_index_ = efind(~sum(ignore_Y_tru_nrm_prm_xj__,1)); %<-- only plot Ys for which both are measured. ;
plot(Y_tru_nrm_prm_xj__(1+0,1+tmp_index_),Y_tru_nrm_prm_xj__(1+1,1+tmp_index_),'ko');
plot(X_est_xt__(1+0,:),X_est_xt__(1+1,:),'rx');
end;%for ni=0:n_i-1;
hold off;
xlabel(u_label_measurement_{1+nv0},'Interpreter','none');
ylabel(u_label_measurement_{1+nv1},'Interpreter','none');
%%%%;
sgtitle(tmp_fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',tmp_fname_fig_pre));
print('-djpeg',tmp_fname_fig_jpg);
close(gcf);
%%%%;
end;%if ~exist(tmp_fname_fig_jpg,'file');
%%%%%%%%
end;%if nshuffle==0;
%%%%%%%%;
tmp_t = toc(tmp_t);
if (flag_verbose>0); disp(sprintf(' %% nshuffle %d/%d SDE_nlp_ijXaABYC_update_all_1: %0.6fs',nshuffle,n_shuffle,tmp_t)); end;
end;%if ~tmp_flag_skip;
%%%%%%%%;
end;%for nshuffle=0:n_shuffle-1+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;end;%for nv0=0:n_var-1; for nv1=nv0+1:n_var-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now collecting the information from each of the variable-pairs. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_x = 2;
dir_PAD_infix_0 = sprintf('%s',str_species);
%%%%;
dir_PAD_output_0_mat = sprintf('%s/dir_PAD_%s_mat',dir_base,dir_PAD_infix_0);
dir_PAD_output_0_jpg = sprintf('%s/dir_PAD_%s_jpg',dir_base,dir_PAD_infix_0);
if ~exist(dir_PAD_output_0_mat,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_0_mat)); mkdir(dir_PAD_output_0_mat); end;
if ~exist(dir_PAD_output_0_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_0_jpg)); mkdir(dir_PAD_output_0_jpg); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_fname_collect_pre = sprintf('%s/%s_collect_0',dir_PAD_output_0_mat,dir_PAD_infix_0);
[tmp_flag_collect_skip,tmp_fname_collect_mat] = open_fname_tmp(tmp_fname_collect_pre);
if flag_recollect | ~tmp_flag_collect_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_t_vv__ = zeros(n_var,n_var);
n_t_use_vv__ = zeros(n_var,n_var);
n_j_vv__ = zeros(n_var,n_var);
n_j_use_vv__ = zeros(n_var,n_var);
n_v_vv__ = zeros(n_var,n_var);
n_v_use_vv__ = zeros(n_var,n_var);
a_est_xasvv_____ = zeros(n_x,n_a,1+n_shuffle,n_var,n_var);
A_est_xxsvv_____ = zeros(n_x,n_x,1+n_shuffle,n_var,n_var);
B_est_omega_svv___ = zeros(1+n_shuffle,n_var,n_var);
B_est_l0_svv___ = zeros(1+n_shuffle,n_var,n_var);
B_est_l1_svv___ = zeros(1+n_shuffle,n_var,n_var);
C_est_omega_svv___ = zeros(1+n_shuffle,n_var,n_var);
C_est_l0_svv___ = zeros(1+n_shuffle,n_var,n_var);
C_est_l1_svv___ = zeros(1+n_shuffle,n_var,n_var);
BtBn_est_xxsvv_____ = zeros(n_x,n_x,1+n_shuffle,n_var,n_var);
CtCn_est_xxsvv_____ = zeros(n_x,n_x,1+n_shuffle,n_var,n_var);
%%%%%%%%%%%%%%%%;
for nv0=0:n_var-1; for nv1=nv0+1:n_var-1;
if (flag_verbose); if ( ((nv0==n_var-2) & (nv1==n_var-1)) | mod(nv1+nv0*n_var,16)==0 ); disp(sprintf(' %% nv0 nv1 (%.2d,%.2d)',nv0,nv1)); end; end;
%%%%%%%%%%%%%%%%;
n_i = n_iid;
n_t_i_ = zeros(n_i,1);
t_it__ = cell(n_i,1);
n_j_i_ = zeros(n_i,1);
n_v_i_ = zeros(n_i,1);
index_nt_from_nj_i__ = cell(n_i,1);
ignore_Y_tru_ixj___ = cell(n_i,1);
Y_tru_ixj___ = cell(n_i,1);
for ni=0:n_i-1;
tmp_retain_Y_xt__ = squeeze(data_measure_iva___(1+ni,1+[nv0,nv1],:));
tmp_ignore_Y_xt__ = ~tmp_retain_Y_xt__;
tmp_Y_xt__= squeeze(data_iva___(1+ni,1+[nv0,nv1],:));
tmp_retain_Y_t_ = reshape(sum(tmp_retain_Y_xt__,1),[n_age,1]);
tmp_index_ = efind(tmp_retain_Y_t_> 0);
t_t_ = age_(1+tmp_index_); n_t = numel(t_t_);
n_j = n_t; %<-- no repeated measurements for PAD data. ;
index_nt_from_nj_ = transpose(0:n_t-1); %<-- no repeated measurements for PAD data. ;
ignore_Y_tru_xj__ = tmp_ignore_Y_xt__(:,1+tmp_index_);
Y_tru_xj__ = tmp_Y_xt__(:,1+tmp_index_);
index_nt_from_nj_i__{1+ni} = index_nt_from_nj_;
ignore_Y_tru_ixj___{1+ni} = ignore_Y_tru_xj__;
Y_tru_ixj___{1+ni} = Y_tru_xj__;
n_j_i_(1+ni) = n_j;
n_t_i_(1+ni) = n_t;
t_it__{1+ni} = t_t_;
n_v_i_(1+ni) = sum(~ignore_Y_tru_xj__,'all');
clear tmp_retain_Y_xt__ tmp_ignore_Y_xt__ tmp_Y_xt__ tmp_retain_Y_t_ tmp_index_ t_t_ n_j index_nt_from_nj_ ignore_Y_tru_xj__ Y_tru_xj__ ;
end;%for ni=0:n_i-1;
n_t_vv__(1+nv0,1+nv1) = sum(n_t_i_);
n_j_vv__(1+nv0,1+nv1) = sum(n_j_i_);
n_v_vv__(1+nv0,1+nv1) = sum(n_v_i_);
n_t_vv__(1+nv1,1+nv0) = sum(n_t_i_);
n_j_vv__(1+nv1,1+nv0) = sum(n_j_i_);
n_v_vv__(1+nv1,1+nv0) = sum(n_v_i_);
%%%%;
% remove monkeys with fewer than 3 measurements. ;
%%%%;
tmp_index_ = efind(n_j_i_>=3);
n_i = numel(tmp_index_);
n_t_i_ = n_t_i_(1+tmp_index_);
t_it__ = t_it__(1+tmp_index_);
n_j_i_ = n_j_i_(1+tmp_index_);
n_v_i_ = n_v_i_(1+tmp_index_);
index_nt_from_nj_i__ = index_nt_from_nj_i__(1+tmp_index_);
ignore_Y_tru_ixj___ = ignore_Y_tru_ixj___(1+tmp_index_);
Y_tru_ixj___ = Y_tru_ixj___(1+tmp_index_);
n_t_use_vv__(1+nv0,1+nv1) = sum(n_t_i_);
n_j_use_vv__(1+nv0,1+nv1) = sum(n_j_i_);
n_v_use_vv__(1+nv0,1+nv1) = sum(n_v_i_);
n_t_use_vv__(1+nv1,1+nv0) = sum(n_t_i_);
n_j_use_vv__(1+nv1,1+nv0) = sum(n_j_i_);
n_v_use_vv__(1+nv1,1+nv0) = sum(n_v_i_);
clear tmp_index_ n_i n_t_i_ t_it__ n_j_i_ n_v_i_ index_nt_from_nj_i__ ignore_Y_tru_ixj___ Y_tru_ixj___;
%%%%;
dir_PAD_infix_1 = sprintf('%s_nv%.2d',str_species,nv0);
dir_PAD_infix_2 = sprintf('%s_nv%.2dnv%.2d',str_species,nv0,nv1);
%%%%;
dir_PAD_output_1_mat = sprintf('%s/dir_PAD_%s_mat',dir_PAD_output_0_mat,dir_PAD_infix_1);
dir_PAD_output_1_jpg = sprintf('%s/dir_PAD_%s_jpg',dir_PAD_output_0_jpg,dir_PAD_infix_1);
if ~exist(dir_PAD_output_1_mat,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_1_mat)); mkdir(dir_PAD_output_1_mat); end;
if ~exist(dir_PAD_output_1_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_1_jpg)); mkdir(dir_PAD_output_1_jpg); end;
%%%%;
dir_PAD_output_2_mat = sprintf('%s/dir_PAD_%s_mat',dir_PAD_output_1_mat,dir_PAD_infix_2);
dir_PAD_output_2_jpg = sprintf('%s/dir_PAD_%s_jpg',dir_PAD_output_1_jpg,dir_PAD_infix_2);
if ~exist(dir_PAD_output_2_mat,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_2_mat)); mkdir(dir_PAD_output_2_mat); end;
if ~exist(dir_PAD_output_2_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_PAD_output_2_jpg)); mkdir(dir_PAD_output_2_jpg); end;
%%%%;
a_est_xas___ = zeros(n_x,n_a,1+n_shuffle);
A_est_xxs___ = zeros(n_x,n_x,1+n_shuffle);
B_est_omega_s_ = zeros(1+n_shuffle,1);
B_est_l0_s_ = zeros(1+n_shuffle,1);
B_est_l1_s_ = zeros(1+n_shuffle,1);
C_est_omega_s_ = zeros(1+n_shuffle,1);
C_est_l0_s_ = zeros(1+n_shuffle,1);
C_est_l1_s_ = zeros(1+n_shuffle,1);
BtBn_est_xxs___ = zeros(n_x,n_x,1+n_shuffle);
CtCn_est_xxs___ = zeros(n_x,n_x,1+n_shuffle);
for nshuffle=0:n_shuffle-1+1;
str_nshuffle = sprintf('s%.4d',nshuffle);
tmp_fname_pre = sprintf('%s/%s_%s',dir_PAD_output_2_mat,dir_PAD_infix_2,str_nshuffle);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre);
if ~exist(tmp_fname_mat,'file'); if (flag_verbose>0); disp(sprintf(' %% %s not found, necessary to create',tmp_fname_mat)); end; end;
if  exist(tmp_fname_mat,'file'); 
if (flag_verbose>1); disp(sprintf(' %% %s found, not necessary to create',tmp_fname_mat)); end; 
tmp_ = load(tmp_fname_mat);
a_est_xas___(:,:,1+nshuffle) = tmp_.a_est_xa__;
A_est_xxs___(:,:,1+nshuffle) = tmp_.A_est_xx__;
B_est_omega_s_(1+nshuffle) = tmp_.B_est_omega;
B_est_l0_s_(1+nshuffle) = tmp_.B_est_l0;
B_est_l1_s_(1+nshuffle) = tmp_.B_est_l1;
C_est_omega_s_(1+nshuffle) = tmp_.C_est_omega;
C_est_l0_s_(1+nshuffle) = tmp_.C_est_l0;
C_est_l1_s_(1+nshuffle) = tmp_.C_est_l1;
BtBn_est_xxs___(:,:,1+nshuffle) = tmp_.BtBn_est_xx__;
CtCn_est_xxs___(:,:,1+nshuffle) = tmp_.CtCn_est_xx__;
clear tmp_;
end;%if  exist(tmp_fname_mat,'file'); 
end;%for nshuffle=0:n_shuffle-1+1;
%%%%;
a_est_xasvv_____(:,:,:,1+nv0,1+nv1) = a_est_xas___;
A_est_xxsvv_____(:,:,:,1+nv0,1+nv1) = A_est_xxs___;
B_est_omega_svv___(:,1+nv0,1+nv1) = B_est_omega_s_;
B_est_l0_svv___(:,1+nv0,1+nv1) = B_est_l0_s_;
B_est_l1_svv___(:,1+nv0,1+nv1) = B_est_l1_s_;
C_est_omega_svv___(:,1+nv0,1+nv1) = C_est_omega_s_;
C_est_l0_svv___(:,1+nv0,1+nv1) = C_est_l0_s_;
C_est_l1_svv___(:,1+nv0,1+nv1) = C_est_l1_s_;
BtBn_est_xxsvv_____(:,:,:,1+nv0,1+nv1) = BtBn_est_xxs___;
CtCn_est_xxsvv_____(:,:,:,1+nv0,1+nv1) = CtCn_est_xxs___;
%%%%;
a_est_xasvv_____(:,:,:,1+nv1,1+nv0) = flip(a_est_xas___,1);
A_est_xxsvv_____(:,:,:,1+nv1,1+nv0) = flip(flip(A_est_xxs___,1),2);
B_est_omega_svv___(:,1+nv1,1+nv0) = B_est_omega_s_;
B_est_l0_svv___(:,1+nv1,1+nv0) = B_est_l0_s_;
B_est_l1_svv___(:,1+nv1,1+nv0) = B_est_l1_s_;
C_est_omega_svv___(:,1+nv1,1+nv0) = C_est_omega_s_;
C_est_l0_svv___(:,1+nv1,1+nv0) = C_est_l0_s_;
C_est_l1_svv___(:,1+nv1,1+nv0) = C_est_l1_s_;
BtBn_est_xxsvv_____(:,:,:,1+nv1,1+nv0) = flip(flip(BtBn_est_xxs___,1),2);
CtCn_est_xxsvv_____(:,:,:,1+nv1,1+nv0) = flip(flip(CtCn_est_xxs___,1),2);
%%%%;
clear a_est_xas___;
clear A_est_xxs___;
clear B_est_omega_s_;
clear B_est_l0_s_;
clear B_est_l1_s_;
clear C_est_omega_s_;
clear C_est_l0_s_;
clear C_est_l1_s_;
clear BtBn_est_xxs___;
clear CtCn_est_xxs___;
%%%%;
%%%%%%%%%%%%%%%%;
end;end;%for nv0=0:n_var-1; for nv1=nv0+1:n_var-1;
%%%%%%%%%%%%%%%%;
save(tmp_fname_collect_mat ...
     ,'n_var','n_u_label_measurement','u_label_measurement_' ...
     ,'str_species' ...
     ,'n_shuffle' ...
     ,'dir_PAD_infix_0','dir_PAD_output_0_mat','dir_PAD_output_0_jpg' ...
     ,'tolerance_master','n_a','flag_regularize_eccentricity_simultaneous','n_iteration_BtBn','MaxFunEvals_use_BtBn','MaxFunEvals_use_simultaneous' ...
     ,'n_t_vv__','n_t_use_vv__' ...
     ,'n_j_vv__','n_j_use_vv__' ...
     ,'n_v_vv__','n_v_use_vv__' ...
     ,'n_x','n_a' ...
     ,'a_est_xasvv_____' ...
     ,'A_est_xxsvv_____' ...
     ,'B_est_omega_svv___' ...
     ,'B_est_l0_svv___' ...
     ,'B_est_l1_svv___' ...
     ,'C_est_omega_svv___' ...
     ,'C_est_l0_svv___' ...
     ,'C_est_l1_svv___' ...
     ,'BtBn_est_xxsvv_____' ...
     ,'CtCn_est_xxsvv_____' ...
     );
close_fname_tmp(tmp_fname_collect_pre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_recollect | ~tmp_flag_collect_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;


