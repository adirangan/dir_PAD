function [output_xyz___,output_l_,x_lim_,y_lim_,z_lim_] = hist3d_0(x_,y_,z_,n_x_bin,n_y_bin,n_z_bin,x_lim_,y_lim_,z_lim_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% function output_l_ = hist3d_0(x_,y_,z_,n_x_bin,n_y_bin,n_z_bin,x_lim_,y_lim_,z_lim_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% 
% Sets up a 3d histogram;
%
% Inputs: 
% x_ double array of length N storing x-coordinates of data. ;
% y_ double array of length N storing y-coordinates of data. ;
% z_ double array of length N storing z-coordinates of data. ;
% x_lim_ double array of length 2 storing min and max values for x-coordinates. ;
% y_lim_ double array of length 2 storing min and max values for y-coordinates. ;
% z_lim_ double array of length 2 storing min and max values for z-coordinates. ;
% n_x_bin integer number of bins in x-direction. ;
% n_y_bin integer number of bins in y-direction. ;
% n_z_bin integer number of bins in z-direction. ;
%
% Outputs:
% output_l_ matrix of size [n_x_bin*n_y_bin*n_z_bin,1] storing 3d-histogram of data. ;
% Note that output_l_ is in sparse format. ;
% output_xyz___ matrix of size [n_x_bin,n_y_bin,n_z_bin] storing 3d-histogram of data. ;
% Note that output_xyz___ is in full format. ;
%
% Test with: hist3d_0();
% 

if (nargin<1);
G1_ = 0.15*randn(1024*4,3) + repmat([-0.15,-0.25,-0.35],1024*4,1);
G2_ = 0.15*randn(1024*2,3) + repmat([+0.15,+0.25,+0.35],1024*2,1);
G3_ = 0.15*randn(1024*1,3) + repmat([+0.45,-0.45,+0.00],1024*1,1);
G_ = [G1_;G2_;G3_];
n_x_bin = 32; n_y_bin = 33; n_z_bin = 34;
h3d_xyz___ = hist3d_0(G_(:,1),G_(:,2),G_(:,3),n_x_bin,n_y_bin,n_z_bin,[-1,+1],[-1,+1],[-1,+1]);
h2d_xy__ = squeeze(sum(h3d_xyz___,3)); h2d_xz__ = squeeze(sum(h3d_xyz___,2)); h2d_yz__ = squeeze(sum(h3d_xyz___,1));
figure(1);clf;figmed;
subplot(1,3,1); imagesc(transpose(h2d_xy__));
axis image; set(gca,'ydir','normal');
set(gca,'XTick',1:n_x_bin); set(gca,'YTick',1:n_y_bin); xtickangle(90);
xlabel('x'); ylabel('y'); title('h2d_xy__','Interpreter','none');
subplot(1,3,2); imagesc(transpose(h2d_xz__));
axis image; set(gca,'ydir','normal'); axisnotick;
set(gca,'XTick',1:n_x_bin); set(gca,'YTick',1:n_z_bin); xtickangle(90);
xlabel('x'); ylabel('z'); title('h2d_xz__','Interpreter','none');
subplot(1,3,3); imagesc(transpose(h2d_yz__));
axis image; set(gca,'ydir','normal'); axisnotick;
set(gca,'XTick',1:n_y_bin); set(gca,'YTick',1:n_z_bin); xtickangle(90);
xlabel('y'); ylabel('z'); title('h2d_yz__','Interpreter','none');
disp('returning');return;
end;%if (nargin<1);

if (nargin<4);
n_x_bin = 32; n_y_bin = 32; nzbin = 33;
end;%if (nargin<4);

if (nargin<7);
x_lim_=[min(x_),max(x_)]; x_lim_ = mean(x_lim_) + diff(x_lim_)/2*1.0625*[-1,1];
y_lim_=[min(y_),max(y_)]; y_lim_ = mean(y_lim_) + diff(y_lim_)/2*1.0625*[-1,1];
z_lim_=[min(z_),max(z_)]; z_lim_ = mean(z_lim_) + diff(z_lim_)/2*1.0625*[-1,1];
end;%if (nargin<7);

x_ = x_(:);
y_ = y_(:);
z_ = z_(:);
n_l = numel(x_);
index_x_l_ = min(n_x_bin-1,max(0,floor(n_x_bin*(x_-min(x_lim_))/(max(x_lim_)-min(x_lim_)))));
index_y_l_ = min(n_y_bin-1,max(0,floor(n_y_bin*(y_-min(y_lim_))/(max(y_lim_)-min(y_lim_)))));
index_z_l_ = min(n_z_bin-1,max(0,floor(n_z_bin*(z_-min(z_lim_))/(max(z_lim_)-min(z_lim_)))));
index_l_ = index_x_l_ + index_y_l_*n_x_bin + index_z_l_*n_x_bin*n_y_bin ;
output_l_ = sparse(1+index_l_,ones(n_l,1),ones(n_l,1),n_x_bin*n_y_bin*n_z_bin,1);
output_xyz___ = reshape(full(output_l_),[n_x_bin,n_y_bin,n_z_bin]);
