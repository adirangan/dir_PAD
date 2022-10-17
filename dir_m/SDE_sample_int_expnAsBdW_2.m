function ...
output_x_ = ...
SDE_sample_int_expnAsBdW_2( ...
 t ...
,Psi_A_xx__ ...
,Lambda_A_x_ ...
,Psi_A_inv_xx__ ...
,B_inv_xx__ ...
,tolerance_SDE ...
);
% samples the function: ;
% \int_{s=0}^{s=t} exp(-A_xx__ * s) * B_inv_xx__ * dW(s) ;
% = \int_{s=0}^{s=t} Psi_A_xx__ * exp(-diag(Lambda_A_x_) * s) * Psi_A_inv_xx__ * B_inv_xx__ * dW(s) ;
% where Psi_A_xx__ * diag(Lambda_A_x_) * Psi_A_inv_xx__ is an eigendecomposition of A_xx__. ;
% Discretizes the total time t into increments dt that are sufficiently small that: ;
%  | exp(-Lambda_A_x_*dt) - (1-Lambda_A_x_*dt) | ;
% is (uniformly) less than tolerance_SDE. ;
% This is done simply by setting dt = sqrt(2*tolerance_SDE)/max(abs(Lambda_A_x_)). ;
% Consequently: ;
% |exp(-Lambda_A_x_ * dt) - (1-Lambda_A_x_*dt)| ;
% \approx | 1 - Lambda_A_x_*dt + 0.5 * (Lambda_A_x_*dt)^.2 - (1-Lambda_A_x_*dt) | ;
% == 0.5 * (Lambda_A_x_*dt).^2 ;
% <= 0.5 * (max(|Lambda_A_x_|)*dt).^2 ;
% == 0.5 * 2 * tolerance_SDE ;

na=0;
if (nargin<1+na); t=[]; end; na=na+1;
if (nargin<1+na); Psi_A_xx__=[]; end; na=na+1;
if (nargin<1+na); Lambda_A_x_=[]; end; na=na+1;
if (nargin<1+na); Psi_A_inv_xx__=[]; end; na=na+1;
if (nargin<1+na); B_inv_xx__=[]; end; na=na+1;
if (nargin<1+na); tolerance_SDE=[]; end; na=na+1;

if isempty(t); t=1; end;
if isempty(Psi_A_xx__); Psi_A_xx__ = eye(2); end;
if isempty(Lambda_A_x_); Lambda_A_x_ = ones(2,1); end;
if isempty(Psi_A_inv_xx__); Psi_A_inv_xx__ = eye(2); end;
if isempty(B_inv_xx__); B_inv_xx__ = eye(2); end;
if isempty(tolerance_SDE); tolerance_SDE = 1e-2; end;

n_var = numel(Lambda_A_x_);

Lambda_A_max = max(abs(Lambda_A_x_));
if (Lambda_A_max<=0); dt_0in = t; end;
if (Lambda_A_max> 0); dt_0in = sqrt(2*tolerance_SDE)/Lambda_A_max; end;
n_dt = max(2,ceil(t/dt_0in)); %<-- take at least 2 steps. ;
dt = t/n_dt;

S_xx__ = Psi_A_inv_xx__*B_inv_xx__ ;
E_xx__ = diag(exp(-Lambda_A_x_*dt));

output_x_ = zeros(n_var,1);
for ndt=0:n_dt-1;
output_x_ = output_x_ + S_xx__ * randn(n_var,1) * sqrt(dt);
S_xx__ = E_xx__ * S_xx__;
end;%for ndt=0:n_dt-1;
output_x_ = Psi_A_xx__ * output_x_;


