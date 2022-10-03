function ddg_ = ddgammaln(x_);
% second derivative of gammaln. ;
% assumes that x_ is real and positive (preferably an integer). ;
if nargin<1;
x_ = [0:0.125:10];
dx = 1e-4;
dg_ = dgammaln(1+x_);
ddg0_ = (gammaln(1+x_+dx) - 2*gammaln(1+x_) + gammaln(1+x_-dx))/dx^2;
ddg1_ = (dgammaln(1+x_+dx) - dgammaln(1+x_-dx))/(2*dx);
ddg2_ = ddgammaln(1+x_);
subplot(1,2,1);
hold on;
plot(1+x_,ddg0_,'ro-');
plot(1+x_,ddg1_,'mo-');
plot(1+x_,ddg2_,'bx-');
hold off;
xlabel('1+x'); ylabel('d^2(gammaln)/dx^2'); title('second derivative of gammaln(1+x)');
legend({'2diff','1diff','approx'},'Location','NorthWest');
subplot(1,2,2);
hold on;
plot(1+x_,ddg0_-ddg2_,'k.-');
hold off;
xlabel('1+x'); ylabel('error'); title('error in second derivative');
figbig;
disp('returning');return;
end;%if nargin<1;

ddg_ = zeros(size(x_));
n_max = ceil(max(x_,[],'all'));
cs_ = cumsum([0,[1:n_max].^-1]) - euler_mascheroni();
dcs_ = diff(cs_);
ddg_(:) = interp1(0.5 + [0:n_max-1],dcs_,x_(:)-1,'spline');
