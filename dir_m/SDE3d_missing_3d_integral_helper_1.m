function ...
[ ...
 Z_base_100_ ...
,l2_stretch_100 ...
,Z_base_110__ ...
,l2_stretch_110 ...
,Z_base_010_ ...
,l2_stretch_010 ...
,Z_base_011__ ...
,l2_stretch_011 ...
,Z_base_001_ ...
,l2_stretch_001 ...
,Z_base_101__ ...
,l2_stretch_101 ...
] = ...
SDE3d_missing_3d_integral_helper_1( ...
 C_omega_ ...
,C_l0 ...
,C_l1 ...
,C_l2 ...
);
%%%%%%%%;
% see test_missing_3d_integral_2.m for an example. ;
%%%%%%%%;

na=0;
if (nargin<1+na); C_omega_=[]; end; na=na+1;
if (nargin<1+na); C_l0=[]; end; na=na+1;
if (nargin<1+na); C_l1=[]; end; na=na+1;
if (nargin<1+na); C_l2=[]; end; na=na+1;
if isempty(C_omega_); C_omega_ = zeros(3,1); end;
if isempty(C_l0); C_l0 = 0; end;
if isempty(C_l1); C_l1 = 0; end;
if isempty(C_l2); C_l2 = 0; end;

%%%%%%%%;
[ ...
 ~ ...
,CtCn_xx__ ...
] = ...
SDE3d_BtBn_1( ...
 [] ...
,C_omega_ ...
,C_l0 ...
,C_l1 ...
,C_l2 ...
);
%%%%%%%%;
[U__,S__,V__] = svd(CtCn_xx__);
S_ = diag(S__);
sigma_ = sqrt(1./max(1e-12,S_));
dsigma__ = diag(sigma_);
isigma__ = inv(dsigma__); %<-- same as below. ;
%%%%%%%%;

%%%%%%%%;
% 100 ; %<-- i.e., only first index given. ;
%%%%%%%%;
%%%%;
W_base_ = transpose(U__)*[1;0;0];
W_perp__ = transpose(U__)*[0 , 0 ; 1 , 0 ; 0 , 1];
SW_base_ = isigma__*W_base_;
SW_perp__ = isigma__*W_perp__;
R_W_perp__ = qr(W_perp__);
R_W_perp__ = R_W_perp__(1:2,1:2);
det_R_W_perp__ = det(R_W_perp__);
[Q_SW_perp__,R_SW_perp__] = qr(SW_perp__); 
Q_SW_perp__ = Q_SW_perp__(:,1:2);
R_SW_perp__ = R_SW_perp__(1:2,1:2);
det_R_SW_perp__ = det(R_SW_perp__);
l2_stretch = (abs(det_R_SW_perp__)/max(1e-12,abs(det_R_W_perp__))).^2;
%disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base_ = SW_base_ - Q_SW_perp__*transpose(Q_SW_perp__)*SW_base_;
l2_stretch_100 = l2_stretch;
Z_base_100_ = Z_base_;
%%%%%%%%;
% 110 ; %<-- i.e., first and second given. ;
%%%%%%%%;
W_base__ = transpose(U__)*[1 , 0 ; 0 , 1 ; 0 , 0];
W_perp_ = transpose(U__)*[0;0;1];
SW_base__ = isigma__*W_base__;
SW_perp_ = isigma__*W_perp_; Q_SW_perp_ = orth(SW_perp_);
l2_stretch = (fnorm(SW_perp_)/max(1e-12,fnorm(W_perp_))).^2;
%disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base__ = SW_base__ - Q_SW_perp_*transpose(Q_SW_perp_)*SW_base__;
l2_stretch_110 = l2_stretch;
Z_base_110__ = Z_base__;
%%%%%%%%;

%%%%%%%%;
% 010 ; %<-- i.e., only second index given. ;
%%%%%%%%;
%%%%;
W_base_ = transpose(U__)*[0;1;0];
W_perp__ = transpose(U__)*[1 , 0 ; 0 , 0 ; 0 , 1];
SW_base_ = isigma__*W_base_;
SW_perp__ = isigma__*W_perp__;
R_W_perp__ = qr(W_perp__);
R_W_perp__ = R_W_perp__(1:2,1:2);
det_R_W_perp__ = det(R_W_perp__);
[Q_SW_perp__,R_SW_perp__] = qr(SW_perp__); 
Q_SW_perp__ = Q_SW_perp__(:,1:2);
R_SW_perp__ = R_SW_perp__(1:2,1:2);
det_R_SW_perp__ = det(R_SW_perp__);
l2_stretch = (abs(det_R_SW_perp__)/max(1e-12,abs(det_R_W_perp__))).^2;
%disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base_ = SW_base_ - Q_SW_perp__*transpose(Q_SW_perp__)*SW_base_;
l2_stretch_010 = l2_stretch;
Z_base_010_ = Z_base_;
%%%%%%%%;
% 011 ; %<-- i.e., second and third given. ;
%%%%%%%%;
W_base__ = transpose(U__)*[0 , 0 ; 1 , 0 ; 0 , 1];
W_perp_ = transpose(U__)*[1;0;0];
SW_base__ = isigma__*W_base__;
SW_perp_ = isigma__*W_perp_; Q_SW_perp_ = orth(SW_perp_);
l2_stretch = (fnorm(SW_perp_)/max(1e-12,fnorm(W_perp_))).^2;
%disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base__ = SW_base__ - Q_SW_perp_*transpose(Q_SW_perp_)*SW_base__;
l2_stretch_011 = l2_stretch;
Z_base_011__ = Z_base__;
%%%%%%%%;

%%%%%%%%;
% 001 ; %<-- i.e., only third index given. ;
%%%%%%%%;
%%%%;
W_base_ = transpose(U__)*[0;0;1];
W_perp__ = transpose(U__)*[1 , 0 ; 0 , 1 ; 0 , 0];
SW_base_ = isigma__*W_base_;
SW_perp__ = isigma__*W_perp__;
R_W_perp__ = qr(W_perp__);
R_W_perp__ = R_W_perp__(1:2,1:2);
det_R_W_perp__ = det(R_W_perp__);
[Q_SW_perp__,R_SW_perp__] = qr(SW_perp__); 
Q_SW_perp__ = Q_SW_perp__(:,1:2);
R_SW_perp__ = R_SW_perp__(1:2,1:2);
det_R_SW_perp__ = det(R_SW_perp__);
l2_stretch = (abs(det_R_SW_perp__)/max(1e-12,abs(det_R_W_perp__))).^2;
%disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base_ = SW_base_ - Q_SW_perp__*transpose(Q_SW_perp__)*SW_base_;
l2_stretch_001 = l2_stretch;
Z_base_001_ = Z_base_;
%%%%%%%%;
% 101 ; %<-- i.e., first and third given. ; (note not periodically arranged). ;
%%%%%%%%;
W_base__ = transpose(U__)*[1 , 0 ; 0 , 0 ; 0 , 1];
W_perp_ = transpose(U__)*[0;1;0];
SW_base__ = isigma__*W_base__;
SW_perp_ = isigma__*W_perp_; Q_SW_perp_ = orth(SW_perp_);
l2_stretch = (fnorm(SW_perp_)/max(1e-12,fnorm(W_perp_))).^2;
%disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base__ = SW_base__ - Q_SW_perp_*transpose(Q_SW_perp_)*SW_base__;
l2_stretch_101 = l2_stretch;
Z_base_101__ = Z_base__;
%%%%%%%%;

