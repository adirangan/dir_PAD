%%%%%%%%;
% testing out: ;
% subjects_measurements_quarterly_50th_percentile_Chimpanzee.xlsx ;
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
fname_xls = sprintf('%s/subjects_measurements_quarterly_50th_percentile_Chimpanzee.xlsx',dir_base);
flag_replot = 1; nf=0;

fname_mat = sprintf('%s/subjects_measurements_quarterly_50th_percentile_Chimpanzee.mat',dir_base);
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
data_measure_iva___ = isfinite(data_iva___) & (data_iva___~=0);
data_measure_i_ = squeeze(sum(sum(data_measure_iva___,3),2));
data_measure_v_ = squeeze(sum(sum(data_measure_iva___,3),1));
data_measure_a_ = squeeze(sum(sum(data_measure_iva___,2),1));
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
%%%%%%%%;

age_lim_ = [min(age_),max(age_)];
nv0 = efind(strcmp(u_label_measurement_,'Body Weight'));
nv1 = efind(strcmp(u_label_measurement_,'Glucose'));
%%%%%%%%;
% Note that nv0 increases roughly linearly until around age 10-15 or so. ;
% plot(age_,squeeze(data_iva___(:,1+[nv0],:)),'o'); grid on;
% This means that a longevity study might start afterwards (e.g., age 15+). ;
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

parameter = struct('type','parameter');
parameter.flag_verbose = 1;
n_var = 2;
data_0in_iva___ = data_iva___(:,1+[nv0,nv1],:);
data_0in_measure_iva___ = data_measure_iva___(:,1+[nv0,nv1],:);
n_q=[];
q_vq__=[];
A_vv__=[];
B_omega=[];
B_l0=[];
B_l1=[];
C_omega=[];
C_l0=[];
C_l1=[];
[ ...
 parameter ...
] = ...
PAD_fit_0( ...
 parameter ...
,n_age ...
,age_ ...
,n_iid ...
,2 ...
,data_0in_iva___ ...
,data_0in_measure_iva___ ...
);
