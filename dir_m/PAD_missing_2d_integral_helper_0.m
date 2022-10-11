function [Z2_base,l2_stretch] = PAD_missing_2d_integral_helper_0(index_missing,C_omega,C_l0,C_l1);
%%%%%%%%;
% if index_missing==0; nlp = 0.5*log(2*pi) - 0.5*(C_l0 + C_l1) + 0.5*log(l2_stretch) + 0.5*Z2_base*XY2.^2; end; %<-- XY1 missing. ;
% if index_missing==1; nlp = 0.5*log(2*pi) - 0.5*(C_l0 + C_l1) + 0.5*log(l2_stretch) + 0.5*Z2_base*XY1.^2; end; %<-- XY2 missing. ;
%%%%%%%%;

na=0;
if (nargin<1+na); index_missing=[]; end; na=na+1;
if (nargin<1+na); C_omega=[]; end; na=na+1;
if (nargin<1+na); C_l0=[]; end; na=na+1;
if (nargin<1+na); C_l1=[]; end; na=na+1;
if isempty(index_missing); index_missing = 0; end;
if isempty(C_omega); C_omega = 0; end;
if isempty(C_l0); C_l0 = 0; end;
if isempty(C_l1); C_l1 = 0; end;

%%%%%%%%;
% first index is missing. ;
%%%%%%%%;
if (index_missing==0);
%W_base_ = [ +sin(C_omega) ; +cos(C_omega) ];
%W_base_perp_ = [ -cos(C_omega) ; +sin(C_omega) ];
%SW_base_ = [ +exp(C_l0/2)*sin(C_omega) ; +exp(C_l1/2)*cos(C_omega) ];
%SW_base_perp_ = [ -exp(C_l0/2)*cos(C_omega) ; +exp(C_l1/2)*sin(C_omega) ];
l2_stretch = exp(C_l0)*cos(C_omega).^2 + exp(C_l1)*sin(C_omega).^2;
%disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
tmp_A = -0.5*sin(2*C_omega)*(exp(C_l0) - exp(C_l1));
tmp_B = exp(C_l0) * cos(C_omega).^2 + exp(C_l1) * sin(C_omega).^2;
tmp_z = tmp_A/max(1e-12,tmp_B);
%Z_base_ = SW_base_ - tmp_A*SW_base_perp_ / fnorm(SW_base_perp_).^2;
%Z_base = ...
%[ ...
%+exp(C_l0/2)*sin(C_omega) + tmp_z * exp(C_l0/2)*cos(C_omega) ; ...
%+exp(C_l1/2)*cos(C_omega) - tmp_z * exp(C_l1/2)*sin(C_omega) ; ...
%];
%Z2_base = fnorm(Z_base_).^2;
Z2_base = ...
+ exp(C_l0)*sin(C_omega).^2 + tmp_z.^2*exp(C_l0)*cos(C_omega).^2 ...
+ exp(C_l1)*cos(C_omega).^2 + tmp_z.^2*exp(C_l1)*sin(C_omega).^2 ...
- (exp(C_l1) - exp(C_l0))*sin(2*C_omega)*tmp_z ...
 ;
end;%if (index_missing==0);


%%%%%%%%;
% second index is missing. ;
%%%%%%%%;
if (index_missing==1);
%W_base_ = [ +cos(C_omega) ; -sin(C_omega) ];
%W_base_perp_ = [ +sin(C_omega) ; +cos(C_omega) ];
%SW_base_ = [ +exp(C_l0/2)*cos(C_omega) ; -exp(C_l1/2)*sin(C_omega) ];
%SW_base_perp_ = [ +exp(C_l0/2)*sin(C_omega) ; +exp(C_l1/2)*cos(C_omega) ];
l2_stretch = exp(C_l0)*sin(C_omega).^2 + exp(C_l1)*cos(C_omega).^2;
%disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
tmp_A = +0.5*sin(2*C_omega)*(exp(C_l0) - exp(C_l1));
tmp_B = exp(C_l0) * sin(C_omega).^2 + exp(C_l1) * cos(C_omega).^2;
tmp_z = tmp_A/max(1e-12,tmp_B);
%Z_base_ = SW_base_ - tmp_A*SW_base_perp_ / fnorm(SW_base_perp_).^2;
%Z_base = ...
%[ ...
%+exp(C_l0/2)*cos(C_omega) - tmp_z * exp(C_l0/2)*sin(C_omega) ; ...
%-exp(C_l1/2)*sin(C_omega) - tmp_z * exp(C_l1/2)*cos(C_omega) ; ...
%];
%Z2_base = fnorm(Z_base_).^2;
Z2_base = ...
+ exp(C_l0)*cos(C_omega).^2 + tmp_z.^2*exp(C_l0)*sin(C_omega).^2 ...
+ exp(C_l1)*sin(C_omega).^2 + tmp_z.^2*exp(C_l1)*cos(C_omega).^2 ...
+ (exp(C_l1) - exp(C_l0))*sin(2*C_omega)*tmp_z ...
 ;
end;%if (index_missing==1);
