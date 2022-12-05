function test_missing_3d_integral_2(rseed);

if nargin<1; rseed=[]; end;
if isempty(rseed); rseed=0; end;

nf=0;

rng(rseed);
n_x = 3;
C_omega_ = 2*pi*rand(3,1);
C_l0 = randn();
C_l1 = randn();
C_l2 = randn();
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
[U__,S__,V__] = svd(CtCn_xx__);
S_ = diag(S__);
sigma_ = sqrt(1./max(1e-12,S_));
dsigma__ = diag(sigma_);
isigma__ = inv(dsigma__); %<-- same as below. ;
%%%%;
P = @(XY0,XY1,XY2) 1/sqrt(2*pi).^3 * sqrt(det(CtCn_xx__)) * exp(-0.5*[XY0 , XY1 , XY2]*CtCn_xx__*[XY0;XY1;XY2]);
%%%%%%%%;
n_xy0 = 128*1.0 + 0;
xy_max = 7.5;
XY0_ = linspace(-xy_max,+xy_max,n_xy0); dl0 = mean(diff(XY0_));
n_xy1 = 128*1.0 + 1;
XY1_ = linspace(-xy_max,+xy_max,n_xy1); dl1 = mean(diff(XY1_));
n_xy2 = 128*1.0 + 2;
XY2_ = linspace(-xy_max,+xy_max,n_xy2); dl2 = mean(diff(XY2_));
%%%%%%%%;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Here assume only first is given. ;
%%%%%%%%;
p_A_ = zeros(n_xy0,1);
nlp_A_ = zeros(n_xy0,1);
for nxy0=0:n_xy0-1;
p_A = 0; XY0 = XY0_(1+nxy0);
for nxy1=0:n_xy1-1;
for nxy2=0:n_xy2-1;
p_A = p_A + dl1*dl2*P(XY0,XY1_(1+nxy1),XY2_(1+nxy2));
end;%for nxy2=0:n_xy2-1;
end;%for nxy1=0:n_xy1-1;
p_A_(1+nxy0) = p_A;
nlp_A_(1+nxy0) = -log(p_A);
end;%for nxy0=0:n_xy0-1;
disp(sprintf(' %% sum(p_A_): %0.16f',sum(p_A_)*dl0));
%%%%%%%%;
%{
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
disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base_ = SW_base_ - Q_SW_perp__*transpose(Q_SW_perp__)*SW_base_;
Z2_base = fnorm(Z_base_).^2;
%}
%%%%;
Z_base_ = Z_base_100_;
l2_stretch = l2_stretch_100;
%%%%;
p_B_ = zeros(n_xy0,1);
nlp_B_ = zeros(n_xy0,1);
for nxy0=0:n_xy0-1;
p_B = 0; XY0 = XY0_(1+nxy0);
Z2 = XY0*transpose(Z_base_)*Z_base_*XY0;
p_B = 1/sqrt(2*pi) * det(isigma__)/max(1e-12,sqrt(l2_stretch)) * exp(-Z2/2);
p_B_(1+nxy0) = p_B;
nlp_B_(1+nxy0) = 0.5*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch) + 0.5*Z2;
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
sgtitle('[1,0,0]');
%%%%%%%%;

%%%%%%%%;
% Here assume first and second components are given. ;
%%%%%%%%;
p_A__ = zeros(n_xy0,n_xy1,1);
nlp_A__ = zeros(n_xy0,n_xy1,1);
for nxy0=0:n_xy0-1;
for nxy1=0:n_xy1-1;
p_A = 0; XY0 = XY0_(1+nxy0); XY1 = XY1_(1+nxy1);
for nxy2=0:n_xy2-1;
p_A = p_A + dl2*P(XY0,XY1,XY2_(1+nxy2));
end;%for nxy2=0:n_xy2-1;
p_A__(1+nxy0,1+nxy1) = p_A;
nlp_A__(1+nxy0,1+nxy1) = -log(p_A);
end;%for nxy1=0:n_xy1-1;
end;%for nxy0=0:n_xy0-1;
disp(sprintf(' %% sum(p_A__): %0.16f',sum(p_A__,'all')*dl0*dl1));
%%%%%%%%;
%{
W_base__ = transpose(U__)*[1 , 0 ; 0 , 1 ; 0 , 0];
W_perp_ = transpose(U__)*[0;0;1];
SW_base__ = isigma__*W_base__;
SW_perp_ = isigma__*W_perp_; Q_SW_perp_ = orth(SW_perp_);
l2_stretch = (fnorm(SW_perp_)/max(1e-12,fnorm(W_perp_))).^2;
disp(sprintf(' %% l2_stretch %0.6f',l2_stretch));
Z_base__ = SW_base__ - Q_SW_perp_*transpose(Q_SW_perp_)*SW_base__;
Z2_base = fnorm(Z_base__).^2;
%}
%%%%;
Z_base__ = Z_base_110__;
l2_stretch = l2_stretch_110;
%%%%;
p_B__ = zeros(n_xy0,n_xy1,1);
nlp_B__ = zeros(n_xy0,n_xy1,1);
for nxy0=0:n_xy0-1;
for nxy1=0:n_xy1-1;
p_B = 0; XY0 = XY0_(1+nxy0); XY1 = XY1_(1+nxy1);
Z2 = [XY0,XY1]*transpose(Z_base__)*Z_base__*[XY0;XY1];
p_B = 1/sqrt(2*pi).^2 * det(isigma__)/max(1e-12,sqrt(l2_stretch)) * exp(-Z2/2);
p_B__(1+nxy0,1+nxy1) = p_B;
nlp_B__(1+nxy0,1+nxy1) = 1.0*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch) + 0.5*Z2;
end;%for nxy1=0:n_xy1-1;
end;%for nxy0=0:n_xy0-1;
%%%%%%%%;
disp(sprintf(' %% p_A_ vs p_B__: %0.16f',fnorm(p_A__-p_B__)./fnorm(p_A__)));
disp(sprintf(' %% nlp_A__ vs nlp_B__: %0.16f',fnorm(nlp_A__-nlp_B__)./fnorm(nlp_A__)));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024]); np=0;
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(p_A__); colorbar;
axis image; xlabel('XY0'); ylabel('XY1'); axisnotick;
title('p_A__','Interpreter','none');
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(p_B__); colorbar;
axis image; xlabel('XY0'); ylabel('XY1'); axisnotick;
title('p_B__','Interpreter','none');
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(nlp_A__); colorbar;
axis image; xlabel('XY0'); ylabel('XY1'); axisnotick;
title('nlp_A__','Interpreter','none');
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(nlp_B__); colorbar;
axis image; xlabel('XY0'); ylabel('XY1'); axisnotick;
title('nlp_B__','Interpreter','none');
%%%%;
sgtitle('[1,1,0]');
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Here assume only second is given. ;
%%%%%%%%;
p_A_ = zeros(n_xy1,1);
nlp_A_ = zeros(n_xy1,1);
for nxy1=0:n_xy1-1;
p_A = 0; XY1 = XY1_(1+nxy1);
for nxy0=0:n_xy0-1;
for nxy2=0:n_xy2-1;
p_A = p_A + dl0*dl2*P(XY0_(1+nxy0),XY1,XY2_(1+nxy2));
end;%for nxy2=0:n_xy2-1;
end;%for nxy0=0:n_xy0-1;
p_A_(1+nxy1) = p_A;
nlp_A_(1+nxy1) = -log(p_A);
end;%for nxy1=0:n_xy1-1;
disp(sprintf(' %% sum(p_A_): %0.16f',sum(p_A_)*dl1));
%%%%%%%%;
Z_base_ = Z_base_010_;
l2_stretch = l2_stretch_010;
%%%%;
p_B_ = zeros(n_xy1,1);
nlp_B_ = zeros(n_xy1,1);
for nxy1=0:n_xy1-1;
p_B = 0; XY1 = XY1_(1+nxy1);
Z2 = XY1*transpose(Z_base_)*Z_base_*XY1;
p_B = 1/sqrt(2*pi) * det(isigma__)/max(1e-12,sqrt(l2_stretch)) * exp(-Z2/2);
p_B_(1+nxy1) = p_B;
nlp_B_(1+nxy1) = 0.5*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch) + 0.5*Z2;
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
sgtitle('[1,0,0]');
%%%%%%%%;

%%%%%%%%;
% Here assume second and third components are given. ;
%%%%%%%%;
p_A__ = zeros(n_xy1,n_xy2,1);
nlp_A__ = zeros(n_xy1,n_xy2,1);
for nxy1=0:n_xy1-1;
for nxy2=0:n_xy2-1;
p_A = 0; XY1 = XY1_(1+nxy1); XY2 = XY2_(1+nxy2);
for nxy0=0:n_xy0-1;
p_A = p_A + dl0*P(XY0_(1+nxy0),XY1,XY2);
end;%for nxy0=0:n_xy0-1;
p_A__(1+nxy1,1+nxy2) = p_A;
nlp_A__(1+nxy1,1+nxy2) = -log(p_A);
end;%for nxy2=0:n_xy2-1;
end;%for nxy1=0:n_xy1-1;
disp(sprintf(' %% sum(p_A__): %0.16f',sum(p_A__,'all')*dl1*dl2));
%%%%%%%%;
Z_base__ = Z_base_011__;
l2_stretch = l2_stretch_011;
%%%%;
p_B__ = zeros(n_xy1,n_xy2,1);
nlp_B__ = zeros(n_xy1,n_xy2,1);
for nxy1=0:n_xy1-1;
for nxy2=0:n_xy2-1;
p_B = 0; XY1 = XY1_(1+nxy1); XY2 = XY2_(1+nxy2);
Z2 = [XY1,XY2]*transpose(Z_base__)*Z_base__*[XY1;XY2];
p_B = 1/sqrt(2*pi).^2 * det(isigma__)/max(1e-12,sqrt(l2_stretch)) * exp(-Z2/2);
p_B__(1+nxy1,1+nxy2) = p_B;
nlp_B__(1+nxy1,1+nxy2) = 1.0*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch) + 0.5*Z2;
end;%for nxy2=0:n_xy2-1;
end;%for nxy1=0:n_xy1-1;
%%%%%%%%;
disp(sprintf(' %% p_A_ vs p_B__: %0.16f',fnorm(p_A__-p_B__)./fnorm(p_A__)));
disp(sprintf(' %% nlp_A__ vs nlp_B__: %0.16f',fnorm(nlp_A__-nlp_B__)./fnorm(nlp_A__)));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024]); np=0;
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(p_A__); colorbar;
axis image; xlabel('XY1'); ylabel('XY2'); axisnotick;
title('p_A__','Interpreter','none');
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(p_B__); colorbar;
axis image; xlabel('XY1'); ylabel('XY2'); axisnotick;
title('p_B__','Interpreter','none');
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(nlp_A__); colorbar;
axis image; xlabel('XY1'); ylabel('XY2'); axisnotick;
title('nlp_A__','Interpreter','none');
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(nlp_B__); colorbar;
axis image; xlabel('XY1'); ylabel('XY2'); axisnotick;
title('nlp_B__','Interpreter','none');
%%%%;
sgtitle('[1,1,0]');
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Here assume only third is given. ;
%%%%%%%%;
p_A_ = zeros(n_xy2,1);
nlp_A_ = zeros(n_xy2,1);
for nxy2=0:n_xy2-1;
p_A = 0; XY2 = XY2_(1+nxy2);
for nxy1=0:n_xy1-1;
for nxy0=0:n_xy0-1;
p_A = p_A + dl0*dl1*P(XY0_(1+nxy0),XY1_(1+nxy1),XY2);
end;%for nxy0=0:n_xy0-1;
end;%for nxy1=0:n_xy1-1;
p_A_(1+nxy2) = p_A;
nlp_A_(1+nxy2) = -log(p_A);
end;%for nxy2=0:n_xy2-1;
disp(sprintf(' %% sum(p_A_): %0.16f',sum(p_A_)*dl2));
%%%%%%%%;
Z_base_ = Z_base_001_;
l2_stretch = l2_stretch_001;
%%%%;
p_B_ = zeros(n_xy2,1);
nlp_B_ = zeros(n_xy2,1);
for nxy2=0:n_xy2-1;
p_B = 0; XY2 = XY2_(1+nxy2);
Z2 = XY2*transpose(Z_base_)*Z_base_*XY2;
p_B = 1/sqrt(2*pi) * det(isigma__)/max(1e-12,sqrt(l2_stretch)) * exp(-Z2/2);
p_B_(1+nxy2) = p_B;
nlp_B_(1+nxy2) = 0.5*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch) + 0.5*Z2;
end;%for nxy2=0:n_xy2-1;
%%%%%%%%;
disp(sprintf(' %% p_A_ vs p_B_: %0.16f',fnorm(p_A_-p_B_)./fnorm(p_A_)));
disp(sprintf(' %% nlp_A_ vs nlp_B_: %0.16f',fnorm(nlp_A_-nlp_B_)./fnorm(nlp_A_)));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024]);
%%%%;
subplot(2,2,1);
plot(XY2_,p_A_,'k-',XY2_,p_B_,'ro'); 
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('probability');
%%%%;
subplot(2,2,2);
plot(XY2_,p_A_./p_B_,'bx-');
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('ratio');
%%%%;
subplot(2,2,3);
plot(XY2_,nlp_A_,'k-',XY2_,nlp_B_,'ro'); 
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('nlp');
%%%%;
subplot(2,2,4);
plot(XY2_,nlp_A_./nlp_B_,'bx-');
xlim(xy_max*[-1,+1]); xlabel('XY');
ylabel('ratio');
%%%%;
sgtitle('[1,0,0]');
%%%%%%%%;

%%%%%%%%;
% Here assume first and third components are given. ;
%%%%%%%%;
p_A__ = zeros(n_xy0,n_xy2,1);
nlp_A__ = zeros(n_xy0,n_xy2,1);
for nxy0=0:n_xy0-1;
for nxy2=0:n_xy2-1;
p_A = 0; XY0 = XY0_(1+nxy0); XY2 = XY2_(1+nxy2);
for nxy1=0:n_xy1-1;
p_A = p_A + dl1*P(XY0,XY1_(1+nxy1),XY2);
end;%for nxy1=0:n_xy1-1;
p_A__(1+nxy0,1+nxy2) = p_A;
nlp_A__(1+nxy0,1+nxy2) = -log(p_A);
end;%for nxy2=0:n_xy2-1;
end;%for nxy0=0:n_xy0-1;
disp(sprintf(' %% sum(p_A__): %0.16f',sum(p_A__,'all')*dl0*dl2));
%%%%%%%%;
Z_base__ = Z_base_101__;
l2_stretch = l2_stretch_101;
%%%%;
p_B__ = zeros(n_xy0,n_xy2,1);
nlp_B__ = zeros(n_xy0,n_xy2,1);
for nxy0=0:n_xy0-1;
for nxy2=0:n_xy2-1;
p_B = 0; XY0 = XY0_(1+nxy0); XY2 = XY2_(1+nxy2);
Z2 = [XY0,XY2]*transpose(Z_base__)*Z_base__*[XY0;XY2];
p_B = 1/sqrt(2*pi).^2 * det(isigma__)/max(1e-12,sqrt(l2_stretch)) * exp(-Z2/2);
p_B__(1+nxy0,1+nxy2) = p_B;
nlp_B__(1+nxy0,1+nxy2) = 1.0*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch) + 0.5*Z2;
end;%for nxy2=0:n_xy2-1;
end;%for nxy0=0:n_xy0-1;
%%%%%%%%;
disp(sprintf(' %% p_A_ vs p_B__: %0.16f',fnorm(p_A__-p_B__)./fnorm(p_A__)));
disp(sprintf(' %% nlp_A__ vs nlp_B__: %0.16f',fnorm(nlp_A__-nlp_B__)./fnorm(nlp_A__)));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024]); np=0;
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(p_A__); colorbar;
axis image; xlabel('XY0'); ylabel('XY2'); axisnotick;
title('p_A__','Interpreter','none');
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(p_B__); colorbar;
axis image; xlabel('XY0'); ylabel('XY2'); axisnotick;
title('p_B__','Interpreter','none');
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(nlp_A__); colorbar;
axis image; xlabel('XY0'); ylabel('XY2'); axisnotick;
title('nlp_A__','Interpreter','none');
%%%%;
subplot(2,2,1+np);np=np+1;
imagesc(nlp_B__); colorbar;
axis image; xlabel('XY0'); ylabel('XY2'); axisnotick;
title('nlp_B__','Interpreter','none');
%%%%;
sgtitle('[1,1,0]');
%%%%%%%%;



