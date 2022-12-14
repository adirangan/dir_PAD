function output_ = gumbel_nll(input_,g_,tol);
if (nargin<1);
rng(1); 
N = 1024*8;
A_ = -log(1-rand(N,1).^0.5);
B_ = (A_ - mean(A_))/std(A_,1);
n_g0 = 1+32; n_g1 = 1+32;
g0_ = linspace(-1,+1,n_g0);
g1_ = linspace(+0,+2,n_g1);
[g0__,g1__] = ndgrid(g0_,g1_);
g_nll__ = zeros(n_g0,n_g1);
for ng0=0:n_g0-1;
for ng1=0:n_g1-1;
tmp_g_ = [g0__(1+ng0,1+ng1),g1__(1+ng0,1+ng1)];
g_nll__(1+ng0,1+ng1) = sum(gumbel_nll(B_,tmp_g_));
end;%for ng1=0:n_g1-1;
end;%for ng0=0:n_g0-1;
[~,opt_ij] = min(g_nll__,[],'all','linear');
[g0_opt_ij,g1_opt_ij] = ind2sub([n_g0,n_g1],opt_ij);
g0_opt = g0__(opt_ij); assert(g0__(opt_ij)==g0_(g0_opt_ij));
g1_opt = g1__(opt_ij); assert(g1__(opt_ij)==g1_(g1_opt_ij));
g_opt_ = [g0_opt,g1_opt];
n_h = 32;
h_x_lim_ = mean(B_) + std(B_,1)*3.5*[-1,+3];
h_x_ = linspace(min(h_x_lim_),max(h_x_lim_),n_h);
h_B_ = hist(B_,h_x_);
h_B_ = h_B_/sum(h_B_)/mean(diff(h_x_));
x_ = linspace(min(h_x_lim_),max(h_x_lim_),1024);
rho_ = gumbel_pdf(x_,g_opt_);
figure(1);clf;figbig;
subplot(1,2,1); 
hold on;
stairs(h_x_,h_B_,'k-','LineWidth',2);
plot(x_,rho_,'g-','LineWidth',2);
hold off;
xlabel('B'); ylabel('#');
xlim([min(h_x_lim_),max(h_x_lim_)]);
title('distribution');
subplot(1,2,2);
hold on;
imagesc(g_nll__,1e5*[0,1]);
set(gca,'XTick',1:n_g1,'XTickLabel',g1_); xtickangle(90);
set(gca,'YTick',1:n_g0,'YTickLabel',g0_);
plot(g1_opt_ij,g0_opt_ij,'go','MarkerFaceColor','g','MarkerSize',8);
set(gca,'FontSize',8);
xlim([1,n_g1]);ylim([1,n_g0]);
hold off;
title('nll');
colormap('hot'); colorbar; xlabel('g1'); ylabel('g0');

disp('returning');return;
end;%if (nargin<1);
if nargin<3; tol = 1e-6; end;
if nargin<2; g_ = [0;1]; end;
mu = g_(1+0);
beta = tol + g_(1+1)^2;
z_ = zeros(size(input_));
z_ = (input_ - mu)/beta;
output_ = zeros(size(input_));
output_(:) = z_(:) + exp(-z_(:)) + log(beta);
