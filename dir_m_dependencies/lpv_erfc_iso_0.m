function lpv = lpv_erfc_iso_0(R,d);
% The integral we are interested in is: ;
% p(lP0) = \int_{\omega(C)} rho(Y_) dY_, ;
% where \omega(C) is the region where nlP0>=nlP0_data. ;
% Now making the transformation X_ = Y_./sigma_quad_, we have: ;
% p(lP0) = prod(sigma_quad_)*\int_{diag(sigmainv_quad_)*\omega(C)} rho(X_) dX_, ;
%        = prod(sigma_quad_)*\int_{B(R)} exp(+nl_I_quad) * exp(-F_nlP0_crit) * exp(-0.5 * (transpose(X_) * X_) ) dX ;
%        = prod(sigma_quad_)*exp(+nl_I_quad)*exp(-F_nlP0_crit) * \int_{B(R)}  * exp(-0.5 * (transpose(X_) * X_) ) dX ;
%        = 1/sqrt(2*pi)^n_d * [Surface of sphere in n_d-dimensions] * \int_{R}^{\infty}  * exp(-0.5 * r^2) r^(n_d-1)dr ;
% where B(R) is now the exterior of the sphere of radius R. ;
% Note that [Surface of sphere in n_d-dimensions] = n_d*pi^(n_d/2) / gamma(1 + n_d/2);
%%%%;
% The integral I(R;n_d) = [ \int_{R}^{\infty} exp(-0.5*r^2) r^(n_d-1) dr ] satisfies: ;
% I(R;1) = sqrt(pi/2) * erfc(R/sqrt(2));
% I(R;2) = exp(-0.5*R^2) ;
% I(R;n_d>=3) = exp(-0.5*R^2)*R^(n_d-2) + (n_d-2)*I(R;n_d-2). ;
%%%%;
% So putting this all together, we have: ;
% p(lP0) = [ n_d*pi^(n_d/2) / gamma(1 + n_d/2) ] * [ 1/sqrt(2*pi)^n_d ] * I(R;n_d) ;
%%%%%%%%;
if (d<=0); lpv = -Inf; lP_0 = -Inf; end;
if (d> 0);
if (mod(d,2)==0); I_init = exp(-0.5*R.^2); tmp_d = 2; end;%if (mod(d,2)==0); %<-- dimension is even. ;
if (mod(d,2)==1); I_init = sqrt(pi/2) * erfc(R./sqrt(2)); tmp_d = 1; end;%if (mod(d,2)==0); %<-- dimension is odd. ;
tmp_I = I_init;
while (tmp_d<d);
tmp_I = tmp_I*tmp_d + exp(-0.5*R.^2).*R.^(tmp_d);
tmp_d = tmp_d+2;
end;%while (tmp_d<d);
lS = log(d) + (d/2)*log(pi) - gammaln(1 + d/2);
lZ = - (d/2)*(log(2) + log(pi));
assert(tmp_d==d);
lpv = log(tmp_I) + lZ + lS;
end;%if (d> 0);

