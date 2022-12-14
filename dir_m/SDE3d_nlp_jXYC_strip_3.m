function ...
[ ...
 nlp_jXYC ...
,CtCn_xx__ ...
] = ...
SDE3d_nlp_jXYC_strip_3( ...
 n_x ...
,n_t ...
,X_xt__ ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,C_omega_ ...
,C_l0 ...
,C_l1 ...
,C_l2 ...
);

if nargin<1;
rng(0);
n_x = 3; n_t = 1024*8;
X_xt__ = randn(n_x,n_t);
n_j = ceil(n_t*2);
index_nt_from_nj_ = max(0,min(n_t-1,floor(n_t*rand(n_j,1))));
ignore_Y_xj__ = rand(n_x,n_j)<=0.25;
Y_xj__ = randn(n_x,n_j);
C_omega_ = 2*pi*rand(3,1);
C_l0 = randn();
C_l1 = randn();
C_l2 = randn();
tmp_t = tic();
[ ...
 nlp_jXYC ...
,CtCn_xx__ ...
] = ...
SDE3d_nlp_jXYC_strip_3( ...
 n_x ...
,n_t ...
,X_xt__ ...
,n_j ...
,index_nt_from_nj_ ...
,ignore_Y_xj__ ...
,Y_xj__ ...
,C_omega_ ...
,C_l0 ...
,C_l1 ...
,C_l2 ...
);
disp(sprintf(' %% nlp_jXYC: %0.6f',nlp_jXYC));
tmp_t = toc(tmp_t); disp(sprintf(' %% SDE3d_nlp_jXYC_strip_3: %0.6fs',tmp_t));
disp(sprintf(' %% returning')); return;
end;%if nargin<1;

flag_verbose=0;

na=0;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); n_t=[]; end; na=na+1;
if (nargin<1+na); X_xt__=[]; end; na=na+1;
if (nargin<1+na); n_j=[]; end; na=na+1;
if (nargin<1+na); index_nt_from_nj_=[]; end; na=na+1;
if (nargin<1+na); ignore_Y_xj__=[]; end; na=na+1;
if (nargin<1+na); Y_xj__=[]; end; na=na+1;
if (nargin<1+na); C_omega_=[]; end; na=na+1;
if (nargin<1+na); C_l0=[]; end; na=na+1;
if (nargin<1+na); C_l1=[]; end; na=na+1;
if (nargin<1+na); C_l2=[]; end; na=na+1;

XY_xj__ = Y_xj__ - X_xt__(:,1+index_nt_from_nj_);

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
if isempty(ignore_Y_xj__); ignore_Y_xj__ = ~isfinite(Y_xj__); end;%if isempty(ignore_Y_xj__); 
XY_xj__(1+efind(ignore_Y_xj__)) = 0;
attend_Y_xj__ = ~ignore_Y_xj__;
attend_Y_111_j_ = transpose((transpose(attend_Y_xj__)*[1;1;1] == 3) & (transpose(attend_Y_xj__)*[0;0;0] == 0));
attend_Y_100_j_ = transpose((transpose(attend_Y_xj__)*[1;0;0] == 1) & (transpose(attend_Y_xj__)*[0;1;1] == 0));
attend_Y_110_j_ = transpose((transpose(attend_Y_xj__)*[1;1;0] == 2) & (transpose(attend_Y_xj__)*[0;0;1] == 0));
attend_Y_010_j_ = transpose((transpose(attend_Y_xj__)*[0;1;0] == 1) & (transpose(attend_Y_xj__)*[1;0;1] == 0));
attend_Y_011_j_ = transpose((transpose(attend_Y_xj__)*[0;1;1] == 2) & (transpose(attend_Y_xj__)*[1;0;0] == 0));
attend_Y_001_j_ = transpose((transpose(attend_Y_xj__)*[0;0;1] == 1) & (transpose(attend_Y_xj__)*[1;1;0] == 0));
attend_Y_101_j_ = transpose((transpose(attend_Y_xj__)*[1;0;1] == 2) & (transpose(attend_Y_xj__)*[0;1;0] == 0));
attend_Y_000_j_ = transpose((transpose(attend_Y_xj__)*[0;0;0] == 0) & (transpose(attend_Y_xj__)*[1;1;1] == 0));
n_111 = sum(attend_Y_111_j_);
n_100 = sum(attend_Y_100_j_);
n_110 = sum(attend_Y_110_j_);
n_010 = sum(attend_Y_010_j_);
n_011 = sum(attend_Y_011_j_);
n_001 = sum(attend_Y_001_j_);
n_101 = sum(attend_Y_101_j_);
n_000 = sum(attend_Y_000_j_);
if (flag_verbose); disp(sprintf(' %% n_j %0.3d n_111 %0.3d n_100 %0.3d n_110 %0.3d n_010 %0.3d n_011 %0.3d n_001 %0.3d n_101 %0.3d n_000 %0.3d',n_j,n_111,n_100,n_110,n_010,n_011,n_001,n_101,n_000)); end;
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
d_111 = + 0.5*3*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) ;
d_100 = + 0.5*1*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch_100) ;
d_110 = + 0.5*2*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch_110) ;
d_010 = + 0.5*1*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch_010) ;
d_011 = + 0.5*2*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch_011) ;
d_001 = + 0.5*1*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch_001) ;
d_101 = + 0.5*2*log(2*pi) - 0.5*(C_l0 + C_l1 + C_l2) + 0.5*log(l2_stretch_101) ;
ZtZn_100__ = zeros(n_x,n_x); ZtZn_100__(1+0,1+0) = transpose(Z_base_100_)*Z_base_100_;
ZtZn_110__ = zeros(n_x,n_x); ZtZn_110__(1+[0,1],1+[0,1]) = transpose(Z_base_110__)*Z_base_110__;
ZtZn_010__ = zeros(n_x,n_x); ZtZn_010__(1+1,1+1) = transpose(Z_base_010_)*Z_base_010_;
ZtZn_011__ = zeros(n_x,n_x); ZtZn_011__(1+[1,2],1+[1,2]) = transpose(Z_base_011__)*Z_base_011__;
ZtZn_001__ = zeros(n_x,n_x); ZtZn_001__(1+2,1+2) = transpose(Z_base_001_)*Z_base_001_;
ZtZn_101__ = zeros(n_x,n_x); ZtZn_101__(1+[0,2],1+[0,2]) = transpose(Z_base_101__)*Z_base_101__;
Z_111 = + 0.5*sum( ( CtCn_xx__*(XY_xj__.*repmat(attend_Y_111_j_,[n_x,1]))).*(XY_xj__.*repmat(attend_Y_111_j_,[n_x,1])) , 'all' ) ;
Z_100 = + 0.5*sum( (ZtZn_100__*(XY_xj__.*repmat(attend_Y_100_j_,[n_x,1]))).*(XY_xj__.*repmat(attend_Y_100_j_,[n_x,1])) , 'all' ) ;
Z_110 = + 0.5*sum( (ZtZn_110__*(XY_xj__.*repmat(attend_Y_110_j_,[n_x,1]))).*(XY_xj__.*repmat(attend_Y_110_j_,[n_x,1])) , 'all' ) ;
Z_010 = + 0.5*sum( (ZtZn_010__*(XY_xj__.*repmat(attend_Y_010_j_,[n_x,1]))).*(XY_xj__.*repmat(attend_Y_010_j_,[n_x,1])) , 'all' ) ;
Z_011 = + 0.5*sum( (ZtZn_011__*(XY_xj__.*repmat(attend_Y_011_j_,[n_x,1]))).*(XY_xj__.*repmat(attend_Y_011_j_,[n_x,1])) , 'all' ) ;
Z_001 = + 0.5*sum( (ZtZn_001__*(XY_xj__.*repmat(attend_Y_001_j_,[n_x,1]))).*(XY_xj__.*repmat(attend_Y_001_j_,[n_x,1])) , 'all' ) ;
Z_101 = + 0.5*sum( (ZtZn_101__*(XY_xj__.*repmat(attend_Y_101_j_,[n_x,1]))).*(XY_xj__.*repmat(attend_Y_101_j_,[n_x,1])) , 'all' ) ;
%%%%%%%%;
nlp_jXYC = ...
+ Z_111 + n_111*d_111 ...
+ Z_100 + n_100*d_100 ...
+ Z_110 + n_110*d_110 ...
+ Z_010 + n_010*d_010 ...
+ Z_011 + n_011*d_011 ...
+ Z_001 + n_001*d_001 ...
+ Z_101 + n_101*d_101 ...
;






