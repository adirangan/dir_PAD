function ...
[ ...
 parameter ...
,nlp_dtXaAB ...
,BtBn_xx__ ...
,nlp_dtXaAB_0 ...
,nlp_dtXaAB_a ...
,nlp_dtXaAB_da_xa__ ...
,nlp_dtXaAB_da_xaxa____ ...
,nlp_dtXaAB_A ...
,nlp_dtXaAB_dA_xx__ ...
,nlp_dtXaAB_dA_xxxx____ ...
,a_opt_xa__ ...
,A_opt_xx__ ...
,Q_xt__ ...
,P_xt__ ...
,dX_xdt__ ...
,dQ_xdt__ ...
,dP_xdt__ ...
,ZX_xdt__ ...
,ZQ_xdt__ ...
,ZP_xdt__ ...
] = ...
SDE_nlp_dtXaAB_0( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
,Q_xt__ ...
,P_xt__ ...
,dX_xdt__ ...
,dQ_xdt__ ...
,dP_xdt__ ...
,ZX_xdt__ ...
,ZQ_xdt__ ...
,ZP_xdt__ ...
);

str_thisfunction = 'SDE_nlp_dtXaAB_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
rng(0);
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
n_t = 128; t_t_ = sort(rand(n_t,1));
n_x = 2; X_xt__ = randn(n_x,n_t);
n_a = 3; a_xa__ = randn(n_x,n_a);
A_xx__ = randn(n_x,n_x);
B_omega = pi/7 ; B_l0 = randn(); B_l1 = randn();
[ ...
 parameter ...
,nlp_dtXaAB ...
,BtBn_xx__ ...
,nlp_dtXaAB_0 ...
,nlp_dtXaAB_a ...
,nlp_dtXaAB_da_xa__ ...
,nlp_dtXaAB_da_xaxa____ ...
,nlp_dtXaAB_A ...
,nlp_dtXaAB_dA_xx__ ...
,nlp_dtXaAB_dA_xxxx____ ...
,a_opt_xa__ ...
,A_opt_xx__ ...
] = ...
SDE_nlp_dtXaAB_0( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
%%%%;
n_test = 8;
scan_max = 0.25;
n_scan = 16+1; scan_ = transpose(linspace(-scan_max,+scan_max,n_scan));
%%%%;
qua_da_sl__ = zeros(n_scan,n_test);
nlp_da_sl__ = zeros(n_scan,n_test);
for ntest=0:n_test-1;
tmp_da_xa__ = randn(n_x,n_a)*fnorm(a_opt_xa__);
for nscan=0:n_scan-1;
tmp_a_xa__ = a_opt_xa__ + tmp_da_xa__*scan_(1+nscan);
qua_da_sl__(1+nscan,1+ntest) = nlp_dtXaAB_0 + nlp_dtXaAB_a ...
+ 1*reshape(nlp_dtXaAB_da_xa__,[1,n_x*n_a])*reshape(tmp_a_xa__,[n_x*n_a,1]) ...
+ 0.5*reshape(tmp_a_xa__,[1,n_x*n_a]) * reshape(nlp_dtXaAB_da_xaxa____,[n_x*n_a,n_x*n_a]) * reshape(tmp_a_xa__,[n_x*n_a,1]) ...
;
[ ...
 parameter ...
,nlp_da_sl__(1+nscan,1+ntest) ...
] = ...
SDE_nlp_dtXaAB_0( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
,n_a ...
,tmp_a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
end;%for nscan=0:n_scan-1;
end;%for ntest=0:n_test-1;
%%%%;
qua_dA_sl__ = zeros(n_scan,n_test);
nlp_dA_sl__ = zeros(n_scan,n_test);
for ntest=0:n_test-1;
tmp_dA_xx__ = randn(n_x,n_x)*fnorm(A_opt_xx__);
for nscan=0:n_scan-1;
tmp_A_xx__ = A_opt_xx__ + tmp_dA_xx__*scan_(1+nscan);
qua_dA_sl__(1+nscan,1+ntest) = nlp_dtXaAB_0 + nlp_dtXaAB_A ...
+ 1*reshape(nlp_dtXaAB_dA_xx__,[1,n_x*n_x])*reshape(tmp_A_xx__,[n_x*n_x,1]) ...
+ 0.5*reshape(tmp_A_xx__,[1,n_x*n_x]) * reshape(nlp_dtXaAB_dA_xxxx____,[n_x*n_x,n_x*n_x]) * reshape(tmp_A_xx__,[n_x*n_x,1]) ...
;
[ ...
 parameter ...
,nlp_dA_sl__(1+nscan,1+ntest) ...
] = ...
SDE_nlp_dtXaAB_0( ...
 parameter ...
,n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
,n_a ...
,a_xa__ ...
,tmp_A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
end;%for nscan=0:n_scan-1;
end;%for ntest=0:n_test-1;
%%%%%%%%;
figure(1);clf;figmed;
%%%%;
subplot(1,2,1);
hold on;
plot(scan_,nlp_da_sl__,'bo');
plot(scan_,qua_da_sl__,'gx');
plot(0,nlp_da_sl__(1+(n_scan-1)/2),'go','MarkerFaceColor',[1,0,0]);
hold off;
xlabel('test');ylabel('nlp'); grid on;
title('nlp_da_sl__','Interpreter','none');
%%%%;
subplot(1,2,2);
hold on;
plot(scan_,nlp_dA_sl__,'bo');
plot(scan_,qua_dA_sl__,'gx');
plot(0,nlp_dA_sl__(1+(n_scan-1)/2),'go','MarkerFaceColor',[1,0,0]);
hold off;
xlabel('test');ylabel('nlp'); grid on;
title('nlp_dA_sl__','Interpreter','none');
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_t=[]; end; na=na+1;
if (nargin<1+na); t_t_=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); X_xt__=[]; end; na=na+1;
if (nargin<1+na); n_a=[]; end; na=na+1;
if (nargin<1+na); a_xa__=[]; end; na=na+1;
if (nargin<1+na); A_xx__=[]; end; na=na+1;
if (nargin<1+na); B_omega=[]; end; na=na+1;
if (nargin<1+na); B_l0=[]; end; na=na+1;
if (nargin<1+na); B_l1=[]; end; na=na+1;
if (nargin<1+na); Q_xt__=[]; end; na=na+1;
if (nargin<1+na); P_xt__=[]; end; na=na+1;
if (nargin<1+na); dX_xdt__=[]; end; na=na+1;
if (nargin<1+na); dQ_xdt__=[]; end; na=na+1;
if (nargin<1+na); dP_xdt__=[]; end; na=na+1;
if (nargin<1+na); ZX_xdt__=[]; end; na=na+1;
if (nargin<1+na); ZQ_xdt__=[]; end; na=na+1;
if (nargin<1+na); ZP_xdt__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;

a_opt_xa__=zeros(n_x,n_a);
A_opt_xx__=zeros(n_x,n_x);

%%%%%%%%;
if (n_t<=1);
nlp_dtXaAB = 0;
[ ...
 ~ ...
,BtBn_xx__ ...
] = ...
SDE_BtBn_0( ...
 [] ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
nlp_dtXaAB_0 = 0;
nlp_dtXaAB_a = 0;
nlp_dtXaAB_da_xa__ = zeros(n_x,n_a);
nlp_dtXaAB_da_xaxa____ = zeros(n_x,n_a,n_x,n_a);
nlp_dtXaAB_A = 0;
nlp_dtXaAB_dA_xx__ = zeros(n_x,n_x);
nlp_dtXaAB_dA_xxxx____ = zeros(n_x,n_x,n_x,n_x);
end;%if (n_t<=1);
%%%%%%%%;
if (n_t> 1);
%%%%;
[ ...
 nlp_dtXaAB ...
,BtBn_xx__ ...
,Q_xt__ ...
,P_xt__ ...
,dX_xdt__ ...
,dQ_xdt__ ...
,dP_xdt__ ...
,ZX_xdt__ ...
,ZQ_xdt__ ...
,ZP_xdt__ ...
] = ...
SDE_nlp_dtXaAB_strip_0( ...
 n_t ...
,t_t_ ...
,n_x ...
,X_xt__ ...
,n_a ...
,a_xa__ ...
,A_xx__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);
%%%%;

if nargout> 3;
%%%%%%%%%%%%%%%%;
n_dt = n_t-1; dt_dt_ = reshape(diff(t_t_),[n_dt,1]); ij_dt_ = 1:n_dt; t_dt_ = t_t_(ij_dt_);
%%%%%%%%%%%%%%%%;
t_ta__ = zeros(n_t,n_a);
t_dta__ = zeros(n_dt,n_a);
dt_dta__ = zeros(n_dt,n_a);
for na=0:n_a-1;
t_ta__(:,1+na) = t_t_.^(1+na-1); 
t_dta__(:,1+na) = t_dt_.^(1+na-1); 
end;%for na=0:n_a-1;
dt_dta__ = diff(t_ta__,1,1);
inv_dt_dt_ = 1./dt_dt_; tmp_ij_ = find(abs(dt_dt_)<1e-12); inv_dt_dt_(tmp_ij_) = sign(dt_dt_(tmp_ij_)).*1e12;
tmp_BtBndt_xxdt___ = 0.5*bsxfun(@times,reshape(BtBn_xx__,[n_x,n_x,1]),reshape(inv_dt_dt_,[1,1,n_dt]));
I_xx__ = eye(n_x,n_x);
nlp_dtXaAB_0 = - n_dt*( (B_l0 + B_l1)/2 - 0.5*n_x*log(2*pi) ) + sum( 0.5*n_x*log(dt_dt_) );
%%%%;
if (flag_verbose>0);
nlp_dtXaAB_P = sum( bsxfun(@times,bsxfun(@times,reshape(ZP_xdt__,[n_x,1,n_dt]),tmp_BtBndt_xxdt___),reshape(ZP_xdt__,[1,n_x,n_dt])) , 'all' );
%nlp_dtXaAB_P = ...
%  +sum( bsxfun(@times,bsxfun(@times,reshape(ZX_xdt__,[n_x,1,n_dt]),tmp_BtBndt_xxdt___),reshape(ZX_xdt__,[1,n_x,n_dt])) , 'all' ) ...
%  -sum( bsxfun(@times,bsxfun(@times,reshape(ZQ_xdt__,[n_x,1,n_dt]),tmp_BtBndt_xxdt___),reshape(ZX_xdt__,[1,n_x,n_dt])) , 'all' ) ...
%  -sum( bsxfun(@times,bsxfun(@times,reshape(ZX_xdt__,[n_x,1,n_dt]),tmp_BtBndt_xxdt___),reshape(ZQ_xdt__,[1,n_x,n_dt])) , 'all' ) ...
%  +sum( bsxfun(@times,bsxfun(@times,reshape(ZQ_xdt__,[n_x,1,n_dt]),tmp_BtBndt_xxdt___),reshape(ZQ_xdt__,[1,n_x,n_dt])) , 'all' ) ...
%  ;
end;%if (flag_verbose>0);
%%%%;
nlp_dtXaAB_a = sum( bsxfun(@times,bsxfun(@times,reshape(ZX_xdt__,[n_x,1,n_dt]),tmp_BtBndt_xxdt___),reshape(ZX_xdt__,[1,n_x,n_dt])) , 'all' );
%{
%%%%%%%%;
% slower than the alternative (below). ;
%%%%%%%%;
tmp_t = tic();
nlp_dtXaAB_da_xaxa____ = zeros(n_x,n_a,n_x,n_a);
for nm0=0:n_x-1;
for nq0=0:n_a-1;
for nm1=0:n_x-1;
for nq1=0:n_a-1;
tmp_RHS_0_xxdt___ = ...
  bsxfun(@times, ...
	 bsxfun(@times, ...
		bsxfun(@times, ...
		       bsxfun(@times, ...
			      bsxfun(@times, ...
				     bsxfun(@times, ...
					    bsxfun(@times, ...
						   ones(n_x,n_x,n_dt) , ...
						   reshape(I_xx__(:,1+nm0),[n_x,1,1]) ) , ...
					    reshape(dt_dta__(:,1+nq0),[1,1,n_dt]) ) , ...
				     reshape(ones(n_dt,1),[1,1,n_dt]) ) , ...
			      reshape(tmp_BtBndt_xxdt___,[n_x,n_x,n_dt]) ) , ...
		       reshape(I_xx__(:,1+nm1),[1,n_x,1]) ) , ...
		reshape(dt_dta__(:,1+nq1),[1,1,n_dt]) ) , ...
	 reshape(ones(n_dt,1),[1,1,n_dt]) ) ...
  ;
tmp_RHS_1_xxdt___ = ...
  bsxfun(@times, ...
	 bsxfun(@times, ...
		bsxfun(@times, ...
		       bsxfun(@times, ...
			      bsxfun(@times, ...
				     bsxfun(@times, ...
					    bsxfun(@times, ...
						   ones(n_x,n_x,n_dt) , ...
						   reshape(A_xx__(:,1+nm0),[n_x,1,1]) ) , ...
					    reshape( t_dta__(:,1+nq0),[1,1,n_dt]) ) , ...
				     reshape(dt_dt_,[1,1,n_dt]) ) , ...
			      reshape(tmp_BtBndt_xxdt___,[n_x,n_x,n_dt]) ) , ...
		       reshape(I_xx__(:,1+nm1),[1,n_x,1]) ) , ...
		reshape(dt_dta__(:,1+nq1),[1,1,n_dt]) ) , ...
	 reshape(ones(n_dt,1),[1,1,n_dt]) ) ...
  ;
tmp_RHS_2_xxdt___ = ...
  bsxfun(@times, ...
	 bsxfun(@times, ...
		bsxfun(@times, ...
		       bsxfun(@times, ...
			      bsxfun(@times, ...
				     bsxfun(@times, ...
					    bsxfun(@times, ...
						   ones(n_x,n_x,n_dt) , ...
						   reshape(I_xx__(:,1+nm0),[n_x,1,1]) ) , ...
					    reshape(dt_dta__(:,1+nq0),[1,1,n_dt]) ) , ...
				     reshape(ones(n_dt,1),[1,1,n_dt]) ) , ...
			      reshape(tmp_BtBndt_xxdt___,[n_x,n_x,n_dt]) ) , ...
		       reshape(A_xx__(:,1+nm1),[1,n_x,1]) ) , ...
		reshape( t_dta__(:,1+nq1),[1,1,n_dt]) ) , ...
	 reshape(dt_dt_,[1,1,n_dt]) ) ...
  ;
tmp_RHS_3_xxdt___ = ...
  bsxfun(@times, ...
	 bsxfun(@times, ...
		bsxfun(@times, ...
		       bsxfun(@times, ...
			      bsxfun(@times, ...
				     bsxfun(@times, ...
					    bsxfun(@times, ...
						   ones(n_x,n_x,n_dt) , ...
						   reshape(A_xx__(:,1+nm0),[n_x,1,1]) ) , ...
					    reshape( t_dta__(:,1+nq0),[1,1,n_dt]) ) , ...
				     reshape(dt_dt_,[1,1,n_dt]) ) , ...
			      reshape(tmp_BtBndt_xxdt___,[n_x,n_x,n_dt]) ) , ...
		       reshape(A_xx__(:,1+nm1),[1,n_x,1]) ) , ...
		reshape( t_dta__(:,1+nq1),[1,1,n_dt]) ) , ...
	 reshape(dt_dt_,[1,1,n_dt]) ) ...
  ;
nlp_dtXaAB_da_xaxa____(1+nm0,1+nq0,1+nm1,1+nq1) = ...
+ 2*sum(tmp_RHS_0_xxdt___,'all') ...
- 2*sum(tmp_RHS_1_xxdt___,'all') ...
- 2*sum(tmp_RHS_2_xxdt___,'all') ...
+ 2*sum(tmp_RHS_3_xxdt___,'all') ...
;
end;%for nq1=0:n_a-1;
end;%for nm1=0:n_x-1;
end;%for nq0=0:n_a-1;
end;%for nm0=0:n_x-1;
%%%%;
nlp_dtXaAB_da_xa__ = zeros(n_x,n_a);
for nm=0:n_x-1;
for nq=0:n_a-1;
tmp_LHS_0_xxdt___ = ...
  bsxfun(@times, ...
	 bsxfun(@times, ...
		bsxfun(@times, ...
		       bsxfun(@times, ...
			      bsxfun(@times, ...
				     ones(n_x,n_x,n_dt) , ...
				     reshape(I_xx__(:,1+nm),[n_x,1,1]) ) , ...
			      reshape(dt_dta__(:,1+nq),[1,1,n_dt]) ) , ...
		       reshape(ones(n_dt,1),[1,1,n_dt]) ) , ...
		reshape(tmp_BtBndt_xxdt___,[n_x,n_x,n_dt]) ) , ...
	 reshape(ZX_xdt__,[1,n_x,n_dt]) ) ...
  ;
tmp_LHS_1_xxdt___ = ...
  bsxfun(@times, ...
	 bsxfun(@times, ...
		bsxfun(@times, ...
		       bsxfun(@times, ...
			      bsxfun(@times, ...
				     ones(n_x,n_x,n_dt) , ...
				     reshape(A_xx__(:,1+nm),[n_x,1,1]) ) , ...
			      reshape( t_dta__(:,1+nq),[1,1,n_dt]) ) , ...
		       reshape(dt_dt_,[1,1,n_dt]) ) , ...
		reshape(tmp_BtBndt_xxdt___,[n_x,n_x,n_dt]) ) , ...
	 reshape(ZX_xdt__,[1,n_x,n_dt]) ) ...
  ;
nlp_dtXaAB_da_xa__(1+nm,1+nq) = ...
- 2*sum(tmp_LHS_0_xxdt___,'all') ...
+ 2*sum(tmp_LHS_1_xxdt___,'all') ...
;
end;%for nq=0:n_a-1;
end;%for nm=0:n_x-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% nlp_dtXaAB_da , nlp_dtXaAB_da %0.6fs',tmp_t)); end;
%}

tmp_t = tic();
%%%%;
nlp_dtXaAB_da_xaxa____ = zeros(n_x,n_a,n_x,n_a);
for nm0=0:n_x-1;
for nq0=0:n_a-1;
for nm1=0:n_x-1;
for nq1=0:n_a-1;
for nj0=0:n_dt-1;
for nn0=0:n_x-1;
for nn1=0:n_x-1;
nlp_dtXaAB_da_xaxa____(1+nm0,1+nq0,1+nm1,1+nq1) = nlp_dtXaAB_da_xaxa____(1+nm0,1+nq0,1+nm1,1+nq1) ...
+ 2*I_xx__(1+nn0,1+nm0) * dt_dta__(1+nj0,1+nq0) *             1 * tmp_BtBndt_xxdt___(1+nn0,1+nn1,1+nj0) * I_xx__(1+nn1,1+nm1) * dt_dta__(1+nj0,1+nq1) *             1 ...
- 2*A_xx__(1+nn0,1+nm0) *  t_dta__(1+nj0,1+nq0) * dt_dt_(1+nj0) * tmp_BtBndt_xxdt___(1+nn0,1+nn1,1+nj0) * I_xx__(1+nn1,1+nm1) * dt_dta__(1+nj0,1+nq1) *             1 ...
- 2*I_xx__(1+nn0,1+nm0) * dt_dta__(1+nj0,1+nq0) *             1 * tmp_BtBndt_xxdt___(1+nn0,1+nn1,1+nj0) * A_xx__(1+nn1,1+nm1) *  t_dta__(1+nj0,1+nq1) * dt_dt_(1+nj0) ...
+ 2*A_xx__(1+nn0,1+nm0) *  t_dta__(1+nj0,1+nq0) * dt_dt_(1+nj0) * tmp_BtBndt_xxdt___(1+nn0,1+nn1,1+nj0) * A_xx__(1+nn1,1+nm1) *  t_dta__(1+nj0,1+nq1) * dt_dt_(1+nj0) ...
;
end;%for nn1=0:n_x-1;
end;%for nn0=0:n_x-1;
end;%for nj0=0:n_dt-1;
end;%for nq1=0:n_a-1;
end;%for nm1=0:n_x-1;
end;%for nq0=0:n_a-1;
end;%for nm0=0:n_x-1;
%%%%;
nlp_dtXaAB_da_xa__ = zeros(n_x,n_a);
for nm0=0:n_x-1;
for nq0=0:n_a-1;
for nj0=0:n_dt-1;
for nn0=0:n_x-1;
for nn1=0:n_x-1;
nlp_dtXaAB_da_xa__(1+nm0,1+nq0) = nlp_dtXaAB_da_xa__(1+nm0,1+nq0) ...
- 2*I_xx__(1+nn0,1+nm0) * dt_dta__(1+nj0,1+nq0) *             1 * tmp_BtBndt_xxdt___(1+nn0,1+nn1,1+nj0) * ZX_xdt__(1+nn1,1+nj0) ...
+ 2*A_xx__(1+nn0,1+nm0) *  t_dta__(1+nj0,1+nq0) * dt_dt_(1+nj0) * tmp_BtBndt_xxdt___(1+nn0,1+nn1,1+nj0) * ZX_xdt__(1+nn1,1+nj0) ...
;
end;%for nn1=0:n_x-1;
end;%for nn0=0:n_x-1;
end;%for nj0=0:n_dt-1;
end;%for nq0=0:n_a-1;
end;%for nm0=0:n_x-1;
%%%%;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% nlp_dtXaAB_da , nlp_dtXaAB_da %0.6fs',tmp_t)); end;

%%%%%%%%%%%%%%%%;
nlp_dtXaAB_A = sum( bsxfun(@times,bsxfun(@times,reshape(dP_xdt__,[n_x,1,n_dt]),tmp_BtBndt_xxdt___),reshape(dP_xdt__,[1,n_x,n_dt])) , 'all' );
tmp_t = tic();
nlp_dtXaAB_dA_xxxx____ = zeros(n_x,n_x,n_x,n_x);
for nx0=0:n_x-1;
for nx1=0:n_x-1;
for nx2=0:n_x-1;
for nx3=0:n_x-1;
nlp_dtXaAB_dA_xxxx____(1+nx0,1+nx1,1+nx2,1+nx3) = 2*sum( 0.5 * BtBn_xx__(1+nx0,1+nx2) * ( transpose(dt_dt_) .* P_xt__(1+nx1,ij_dt_) .* P_xt__(1+nx3,ij_dt_) ) );
%nlp_dtXaAB_dA_xxxx____(1+nx0,1+nx1,1+nx2,1+nx3) = 2*sum( reshape(tmp_BtBndt_xxdt___(1+nx0,1+nx2,:),[1,n_dt]) .* reshape(dt_dt_.^2,[1,n_dt]) .* P_xt__(1+nx1,ij_dt_) .* P_xt__(1+nx3,ij_dt_) );
end;%for nx3=0:n_x-1;
end;%for nx2=0:n_x-1;
end;%for nx1=0:n_x-1;
end;%for nx0=0:n_x-1;
nlp_dtXaAB_dA_xx__ = zeros(n_x,n_x);
for nx0=0:n_x-1;
for nx1=0:n_x-1;
for nx2=0:n_x-1;
nlp_dtXaAB_dA_xx__(1+nx0,1+nx1) = nlp_dtXaAB_dA_xx__(1+nx0,1+nx1) - 2*sum( 0.5 * BtBn_xx__(1+nx0,1+nx2) * ( P_xt__(1+nx1,ij_dt_) .* dP_xdt__(1+nx2,ij_dt_) ) );
%nlp_dtXaAB_dA_xx__(1+nx0,1+nx1) = nlp_dtXaAB_dA_xx__(1+nx0,1+nx1) + 2*sum( reshape(tmp_BtBndt_xxdt___(1+nx0,1+nx2,:),[1,n_dt]) .* reshape(dt_dt_,[1,n_dt]) .* P_xt__(1+nx1,ij_dt_) .* dP_xdt__(1+nx2,ij_dt_) );
end;%for nx2=0:n_x-1;
end;%for nx1=0:n_x-1;
end;%for nx0=0:n_x-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% nlp_dtXaAB_dA , nlp_dtXaAB_dA %0.6fs',tmp_t)); end;
%%%%%%%%%%%%%%%%;
if (flag_verbose>0);
tmp_nlp_P = nlp_dtXaAB_0 + nlp_dtXaAB_P;
tmp_nlp_a = nlp_dtXaAB_0 + nlp_dtXaAB_a ...
+ 1*reshape(nlp_dtXaAB_da_xa__,[1,n_x*n_a])*reshape(a_xa__,[n_x*n_a,1]) ...
+ 0.5*reshape(a_xa__,[1,n_x*n_a]) * reshape(nlp_dtXaAB_da_xaxa____,[n_x*n_a,n_x*n_a]) * reshape(a_xa__,[n_x*n_a,1]) ...
;
tmp_nlp_A = nlp_dtXaAB_0 + nlp_dtXaAB_A ...
+ 1*reshape(nlp_dtXaAB_dA_xx__,[1,n_x*n_x])*reshape(A_xx__,[n_x*n_x,1]) ...
+ 0.5*reshape(A_xx__,[1,n_x*n_x]) * reshape(nlp_dtXaAB_dA_xxxx____,[n_x*n_x,n_x*n_x]) * reshape(A_xx__,[n_x*n_x,1]) ...
;
disp(sprintf(' %% nlp_dtXaAB vs tmp_nlp_P: %0.16f',fnorm(nlp_dtXaAB-tmp_nlp_P)/fnorm(nlp_dtXaAB)));
disp(sprintf(' %% nlp_dtXaAB vs tmp_nlp_a: %0.16f',fnorm(nlp_dtXaAB-tmp_nlp_a)/fnorm(nlp_dtXaAB)));
disp(sprintf(' %% nlp_dtXaAB vs tmp_nlp_A: %0.16f',fnorm(nlp_dtXaAB-tmp_nlp_A)/fnorm(nlp_dtXaAB)));

end;%if (flag_verbose>0);
%%%%%%%%%%%%%%%%;
end;%if nargout> 3;

if nargout> 7;
a_opt_xa__ = reshape( - pinv(reshape(nlp_dtXaAB_da_xaxa____,[n_x*n_a,n_x*n_a]),tolerance_master) * reshape(nlp_dtXaAB_da_xa__,[n_x*n_a,1]) , [n_x,n_a] );
A_opt_xx__ = reshape( - pinv(reshape(nlp_dtXaAB_dA_xxxx____,[n_x^2,n_x^2]),tolerance_master) * reshape(nlp_dtXaAB_dA_xx__,[n_x^2,1]) , [n_x,n_x] );
end;%if nargout> 7;

%%%%;
end;%if (n_t> 1);
%%%%%%%%;






