function test_missing_2d_integral_0(rseed);

if nargin<1; rseed=[]; end;
if isempty(rseed); rseed=0; end;

nf=0;

rng(rseed);
n_x = 2;
C_omega = 2*pi*rand();
C_l0 = randn();
C_l1 = randn();
[ ...
 ~ ...
,CtCn_xx__ ...
] = ...
PAD_BtBn_0( ...
 [] ...
,C_omega ...
,C_l0 ...
,C_l1 ...
);
[U__,S__,V__] = svd(CtCn_xx__);
S_ = diag(S__);
sigma_ = sqrt(1./max(1e-12,S_));
dsigma__ = diag(sigma_);
isigma__ = inv(dsigma__); %<-- same as below. ;
%%%%;
U__ = [ +cos(C_omega) , -sin(C_omega) ; +sin(C_omega) , +cos(C_omega) ];
isigma__ = diag(exp([C_l0/2 , C_l1/2]));
%%%%;
P = @(XY0,XY1) 1/(2*pi) * sqrt(det(CtCn_xx__)) * exp(-0.5*[XY0 , XY1]*CtCn_xx__*[XY0;XY1]);
%%%%%%%%;
n_xy0 = 1024*1.0 + 0;
xy_max = 7.5;
XY0_ = linspace(-xy_max,+xy_max,n_xy0); dl0 = mean(diff(XY0_));
n_xy1 = 1024*1.0 + 1;
XY1_ = linspace(-xy_max,+xy_max,n_xy1); dl1 = mean(diff(XY1_));
%%%%%%%%;

%%%%%%%%;
% Here assume second component is missing. ;
%%%%%%%%;
p_A_ = zeros(n_xy0,1);
nlp_A_ = zeros(n_xy0,1);
for nxy0=0:n_xy0-1;
p_A = 0; XY0 = XY0_(1+nxy0);
for nxy1=0:n_xy1-1;
p_A = p_A + dl1*P(XY0,XY1_(1+nxy1));
end;%for nxy1=0:n_xy1-1;
p_A_(1+nxy0) = p_A;
nlp_A_(1+nxy0) = -log(p_A);
end;%for nxy0=0:n_xy0-1;
disp(sprintf(' %% sum(p_A_): %0.16f',sum(p_A_)*dl0));
%%%%%%%%;
W_base_ = transpose(U__)*[1;0];
W_base_perp_ = [-W_base_(2) ; +W_base_(1)];
SW_base_ = isigma__*W_base_;
SW_base_perp_ = isigma__*W_base_perp_;
l2_stretch = (fnorm(SW_base_perp_)/max(1e-12,fnorm(W_base_perp_))).^2;
disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base_ = SW_base_ - dot(SW_base_perp_,SW_base_)*SW_base_perp_ / fnorm(SW_base_perp_).^2;
Z2_base = fnorm(Z_base_).^2;
%%%%;
W_base_ = [ +cos(C_omega) ; -sin(C_omega) ];
W_base_perp_ = [ +sin(C_omega) ; +cos(C_omega) ];
SW_base_ = [ +exp(C_l0/2)*cos(C_omega) ; -exp(C_l1/2)*sin(C_omega) ];
SW_base_perp_ = [ +exp(C_l0/2)*sin(C_omega) ; +exp(C_l1/2)*cos(C_omega) ];
l2_stretch = exp(C_l0)*sin(C_omega).^2 + exp(C_l1)*cos(C_omega).^2;
disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
tmp_A = +0.5*sin(2*C_omega)*(exp(C_l0) - exp(C_l1));
tmp_B = exp(C_l0) * sin(C_omega).^2 + exp(C_l1) * cos(C_omega).^2;
tmp_z = tmp_A/max(1e-12,tmp_B);
%Z_base_ = SW_base_ - tmp_A*SW_base_perp_ / fnorm(SW_base_perp_).^2;
Z_base = ...
[ ...
+exp(C_l0/2)*cos(C_omega) - tmp_z * exp(C_l0/2)*sin(C_omega) ; ...
-exp(C_l1/2)*sin(C_omega) - tmp_z * exp(C_l1/2)*cos(C_omega) ; ...
];
%Z2_base = fnorm(Z_base_).^2;
Z2_base = ...
+ exp(C_l0)*cos(C_omega).^2 + tmp_z.^2*exp(C_l0)*sin(C_omega).^2 ...
+ exp(C_l1)*sin(C_omega).^2 + tmp_z.^2*exp(C_l1)*cos(C_omega).^2 ...
+ (exp(C_l1) - exp(C_l0))*sin(2*C_omega)*tmp_z ...
 ;
[Z2_base,l2_stretch] = PAD_missing_2d_integral_helper_0(1,C_omega,C_l0,C_l1);
%%%%;
p_B_ = zeros(n_xy0,1);
nlp_B_ = zeros(n_xy0,1);
for nxy0=0:n_xy0-1;
p_B = 0; XY0 = XY0_(1+nxy0);
Z2 = Z2_base*XY0.^2;
p_B = 1/sqrt(2*pi) * det(isigma__)/max(1e-12,sqrt(l2_stretch)) * exp(-Z2/2);
p_B_(1+nxy0) = p_B;
nlp_B_(1+nxy0) = 0.5*log(2*pi) - 0.5*(C_l0 + C_l1) + 0.5*log(l2_stretch) + 0.5*Z2_base*XY0.^2;
end;%for nxy0=0:n_xy0-1;
%%%%%%%%;
disp(sprintf(' %% p_A_ vs p_B_: %0.16f',fnorm(p_A_-p_B_)./fnorm(p_A_)));
disp(sprintf(' %% nlp_A_ vs nlp_B_: %0.16f',fnorm(nlp_A_-nlp_B_)./fnorm(nlp_A_)));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024]);
%%%%;
subplot(2,2,1);
plot(XY0_,p_A_,'k-',XY0_,p_B_,'ro'); 
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('probability');
%%%%;
subplot(2,2,2);
plot(XY0_,p_A_./p_B_,'bx-');
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('ratio');
%%%%;
subplot(2,2,3);
plot(XY0_,nlp_A_,'k-',XY0_,nlp_B_,'ro'); 
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('nlp');
%%%%;
subplot(2,2,4);
plot(XY0_,nlp_A_./nlp_B_,'bx-');
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('ratio');
%%%%;
sgtitle('[1,0]');
%%%%%%%%;

%%%%%%%%;
% Here assume first component is missing. ;
%%%%%%%%;
p_A_ = zeros(n_xy1,1);
nlp_A_ = zeros(n_xy1,1);
for nxy1=0:n_xy1-1;
p_A = 0; XY1 = XY1_(1+nxy1);
for nxy0=0:n_xy0-1;
p_A = p_A + dl0*P(XY0_(1+nxy0),XY1);
end;%for nxy0=0:n_xy0-1;
p_A_(1+nxy1) = p_A;
nlp_A_(1+nxy1) = -log(p_A);
end;%for nxy1=0:n_xy1-1;
disp(sprintf(' %% sum(p_A_): %0.16f',sum(p_A_)*dl1));
%%%%%%%%;
W_base_ = transpose(U__)*[0;1];
W_base_perp_ = [-W_base_(2) ; +W_base_(1)];
SW_base_ = isigma__*W_base_;
SW_base_perp_ = isigma__*W_base_perp_;
l2_stretch = (fnorm(SW_base_perp_)/max(1e-12,fnorm(W_base_perp_))).^2;
disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base_ = SW_base_ - dot(SW_base_perp_,SW_base_)*SW_base_perp_ / fnorm(SW_base_perp_).^2;
Z2_base = fnorm(Z_base_).^2;
%%%%;
W_base_ = [ +sin(C_omega) ; +cos(C_omega) ];
W_base_perp_ = [ -cos(C_omega) ; +sin(C_omega) ];
SW_base_ = [ +exp(C_l0/2)*sin(C_omega) ; +exp(C_l1/2)*cos(C_omega) ];
SW_base_perp_ = [ -exp(C_l0/2)*cos(C_omega) ; +exp(C_l1/2)*sin(C_omega) ];
l2_stretch = exp(C_l0)*cos(C_omega).^2 + exp(C_l1)*sin(C_omega).^2;
disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
tmp_A = -0.5*sin(2*C_omega)*(exp(C_l0) - exp(C_l1));
tmp_B = exp(C_l0) * cos(C_omega).^2 + exp(C_l1) * sin(C_omega).^2;
tmp_z = tmp_A/max(1e-12,tmp_B);
%Z_base_ = SW_base_ - tmp_A*SW_base_perp_ / fnorm(SW_base_perp_).^2;
Z_base = ...
[ ...
+exp(C_l0/2)*sin(C_omega) + tmp_z * exp(C_l0/2)*cos(C_omega) ; ...
+exp(C_l1/2)*cos(C_omega) - tmp_z * exp(C_l1/2)*sin(C_omega) ; ...
];
%Z2_base = fnorm(Z_base_).^2;
Z2_base = ...
+ exp(C_l0)*sin(C_omega).^2 + tmp_z.^2*exp(C_l0)*cos(C_omega).^2 ...
+ exp(C_l1)*cos(C_omega).^2 + tmp_z.^2*exp(C_l1)*sin(C_omega).^2 ...
- (exp(C_l1) - exp(C_l0))*sin(2*C_omega)*tmp_z ...
 ;
[Z2_base,l2_stretch] = PAD_missing_2d_integral_helper_0(0,C_omega,C_l0,C_l1);
%%%%;
p_B_ = zeros(n_xy1,1);
nlp_B_ = zeros(n_xy1,1);
for nxy1=0:n_xy1-1;
p_B = 0; XY1 = XY1_(1+nxy1);
Z2 = Z2_base*XY1.^2;
p_B = 1/sqrt(2*pi) * det(isigma__)/max(1e-12,sqrt(l2_stretch)) * exp(-Z2/2);
p_B_(1+nxy1) = p_B;
nlp_B_(1+nxy1) = 0.5*log(2*pi) - 0.5*(C_l0 + C_l1) + 0.5*log(l2_stretch) + 0.5*Z2_base*XY1.^2;
end;%for nxy1=0:n_xy1-1;
%%%%%%%%;
disp(sprintf(' %% p_A_ vs p_B_: %0.16f',fnorm(p_A_-p_B_)./fnorm(p_A_)));
disp(sprintf(' %% nlp_A_ vs nlp_B_: %0.16f',fnorm(nlp_A_-nlp_B_)./fnorm(nlp_A_)));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024]);
%%%%;
subplot(2,2,1);
plot(XY1_,p_A_,'k-',XY1_,p_B_,'ro'); 
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('probability');
%%%%;
subplot(2,2,2);
plot(XY1_,p_A_./p_B_,'bx-');
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('ratio');
%%%%;
subplot(2,2,3);
plot(XY1_,nlp_A_,'k-',XY1_,nlp_B_,'ro'); 
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('nlp');
%%%%;
subplot(2,2,4);
plot(XY1_,nlp_A_./nlp_B_,'bx-');
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('ratio');
%%%%;
sgtitle('[0,1]');
%%%%%%%%;


