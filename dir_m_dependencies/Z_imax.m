function [zone_max,zone_max_ij] = Z_imax(verbose,Z_,Z_min);
% finds interior maximum of Z_;
if nargin<1;
verbose=1;
n_x = 1024; x_ = linspace(0,2*pi,n_x);
figure(1);clf;figbig;
prows=4;pcols=5;np=1;
Z_ = exp(0.00-x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(x_-2*pi); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(0.00-x_).*sin(x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(x_-2*pi).*sin(x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(0.00-x_).*cos(x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(x_-2*pi).*cos(x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = cos(x_)+2; [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = (x_-pi).^2; [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(0.00-x_)+0.15*cos(3*x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(x_-2*pi)+0.15*cos(3*x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(0.00-x_)+0.15*sin(3*x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(x_-2*pi)+0.15*sin(3*x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(0.00-x_)+0.1*sin(3*x_)+1; [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(x_-2*pi)+0.1*sin(3*x_)+1; [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(0.00-x_)+0.1*cos(3*x_)+1; [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(x_-2*pi)+0.1*cos(3*x_)+1; [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(0.00-x_)+0.5*cos(3*x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(x_-2*pi)+0.5*cos(3*x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(0.00-x_)+0.5*sin(3*x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
Z_ = exp(x_-2*pi)+0.5*sin(3*x_); [Z,ij] = Z_imax(0,Z_,0); subplot(prows,pcols,np); plot(x_,Z_,'k.-',x_(ij),Z,'ro',x_,zeros(n_x,1),'k:'); xlim([0,2*pi]); np=np+1; disp(sprintf(" %% Z %0.16f ij-1 %d",Z,ij-1));
figbig;
sgtitle('Z_imax');
disp('returning'); return;
end;%if nargin<1;

if (verbose);
disp(sprintf(' %% [entering Z_imax]'));
if (verbose>1); disp(sprintf(' %% Z_min: %0.16f; Z_: ',Z_min)); disp(num2str(transpose(Z_(:)))); end;
end;%if (verbose);
n_Z = numel(Z_);
[Z_max,ij] = max(Z_); zone_max = Z_max; zone_max_ij = ij;
if (Z_max<Z_min); zone_max = Z_max; zone_max_ij = []; end;
if (Z_max>=Z_min);
ij_sub_ = find(Z_>=Z_min);
n_sub = length(ij_sub_);
n_zone = 1;
for nsub=1:n_sub-1;
if (ij_sub_(nsub+1)-ij_sub_(nsub))>1; n_zone = n_zone+1; end;
end;%for nsub=1:n_sub-1;
if (verbose); disp(sprintf(' %% found %d zones',n_zone)); end;
zone_start_ = ones(1+n_zone,1); zone_start_(1) = ij_sub_(1);
zone_end_ = ones(1+n_zone,1); zone_end_(n_zone) = ij_sub_(end);
nzone = 1;
for nsub=1:n_sub-1;
if (ij_sub_(nsub+1)-ij_sub_(nsub))>1; 
zone_end_(nzone) = ij_sub_(nsub); 
zone_start_(1+nzone) = ij_sub_(nsub+1); 
nzone = nzone+1; 
end;%if (ij_sub_(nsub+1)-ij_sub_(nsub))>1; 
end;%for nsub=1:n_sub-1;
if (verbose); 
for nzone=1:n_zone;
disp(sprintf(' %% zone %d: [%d:%d]',nzone,zone_start_(nzone),zone_end_(nzone)));
end;%for nzone=1:n_zone;
end;%if (verbose); 
zone_ij_ = cell(n_zone,1);
zone_numel_ = zeros(n_zone,1);
for nzone=1:n_zone;
zone_ij_{nzone} = [zone_start_(nzone):zone_end_(nzone)];
zone_numel_(nzone) = numel(zone_ij_{nzone});
end;%for nzone=1:n_zone;
zone_max_ = zeros(n_zone,1);
zone_max_ij_ = zeros(n_zone,1);
zone_valid_ = zeros(n_zone,1);
for nzone=1:n_zone;
[zone_max_(nzone),tmp_ij] = max(Z_(zone_ij_{nzone}));
zone_max_ij_(nzone) = zone_ij_{nzone}(tmp_ij);
if (nzone>1 & nzone<n_zone); zone_valid_(nzone) = 1; end;
if (nzone==1 & zone_max_ij_(nzone)>1 & zone_max_ij_(nzone)<n_Z); zone_valid_(nzone) = 1; end;
if (nzone==n_zone & zone_max_ij_(nzone)>1 & zone_max_ij_(nzone)<n_Z); zone_valid_(nzone) = 1; end;
end;%for nzone=1:n_zone;
if (verbose); 
for nzone=1:n_zone;
disp(sprintf(' %% zone %d: [%d:%d], max %0.2f, ij %d, valid %d',nzone,zone_start_(nzone),zone_end_(nzone),zone_max_(nzone),zone_max_ij_(nzone),zone_valid_(nzone)));
end;%for nzone=1:n_zone;
end;%if (verbose); 
if (sum(zone_valid_));
if (verbose); disp(sprintf(' %% found %d valid zones',sum(zone_valid_))); end;
tmp_valid_ij_ = find(zone_valid_);
[zone_max,z_ij] = max(zone_max_(tmp_valid_ij_));
z_ij = tmp_valid_ij_(z_ij);
if (verbose); disp(sprintf(' %% found max in zone z_ij: %d',z_ij)); end;
assert(zone_max==zone_max_(z_ij));
zone_max_ij = zone_max_ij_(z_ij);
end;%if (sum(zone_valid_));
if (sum(zone_valid_)==0);
if (verbose); disp(sprintf(' %% found no valid zones')); end;
[zone_max,z_ij] = max(zone_max_);
assert(zone_max==zone_max_(z_ij));
zone_max_ij = zone_max_ij_(z_ij);
if (verbose); disp(sprintf(' %% found max in zone z_ij: %d',zone_max_ij)); end;
if (0);
elseif ((zone_max_ij==1 | zone_max_ij==n_Z) & n_zone==1 & zone_numel_(1)==n_Z);
tmp_min = min(Z_) + max(1e-12,abs(min(Z_))*1e-12);
if (verbose); disp(sprintf(' %% only one big zone; rerunning with new minimum Z_min %0.16f --> %0.16f',min(Z_),tmp_min)); end;
[zone_max,zone_max_ij] = Z_imax(verbose,Z_,tmp_min);
elseif (zone_max_ij==1 & zone_numel_(z_ij)<n_Z);
if (verbose); disp(sprintf(' %% zone on left side')); end;
zone_lmax_ij_ = localmax(Z_(zone_ij_{z_ij})); zone_lmax_ij_ = zone_ij_{z_ij}(zone_lmax_ij_);
if numel(zone_lmax_ij_)>1; 
if (verbose); disp(sprintf(' %% found local maximum')); end;
zone_max_ij = zone_lmax_ij_(2); zone_max = Z_(zone_max_ij);
end;%if numel(zone_lmax_ij_)>1; 
if numel(zone_lmax_ij_)<=1; 
if (verbose); disp(sprintf(' %% could not find local maximum, choosing midpoint.')); end;
tmp_z = 0.5*(max(Z_(zone_ij_{z_ij})) + min(Z_(zone_ij_{z_ij})));
[~,zone_max_ij] = min(abs(Z_(zone_ij_{z_ij})-tmp_z));
zone_max_ij = zone_ij_{z_ij}(zone_max_ij);
zone_max = Z_(zone_max_ij);
end;%if numel(zone_lmax_ij_)<=1; 
elseif (zone_max_ij==n_Z & zone_numel_(z_ij)<n_Z);
if (verbose); disp(sprintf(' %% zone on right side')); end;
zone_lmax_ij_ = localmax(Z_(zone_ij_{z_ij})); zone_lmax_ij_ = zone_ij_{z_ij}(zone_lmax_ij_);
if numel(zone_lmax_ij_)>1; 
if (verbose); disp(sprintf(' %% found local maximum')); end;
zone_max_ij = zone_lmax_ij_(2); zone_max = Z_(zone_max_ij);
end;%if numel(zone_lmax_ij_)>1; 
if numel(zone_lmax_ij_)<=1; 
if (verbose); disp(sprintf(' %% could not find local maximum, choosing midpoint.')); end;
tmp_z = 0.5*(max(Z_(zone_ij_{z_ij})) + min(Z_(zone_ij_{z_ij})));
[~,zone_max_ij] = min(abs(Z_(zone_ij_{z_ij})-tmp_z));
zone_max_ij = zone_ij_{z_ij}(zone_max_ij);
zone_max = Z_(zone_max_ij);
end;%if numel(zone_lmax_ij_)<=1; 
end;%if (0);
end;%if (sum(zone_valid_)==0);
if (verbose); disp(sprintf(' %% zone_max %0.2f<--%0.2f zone_max_ij %d',zone_max,Z_(zone_max_ij),zone_max_ij)); end;
end;%if (Z_max>=Z_min);
if (verbose);disp(sprintf(' %% [finished Z_imax]')); end;