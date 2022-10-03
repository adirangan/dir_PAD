function ...
[ ...
 parameter ...
,C_n_xx__ ...
,omega ...
,l0 ...
,l1 ...
,niteration_newton ...
,nlp ...
,dw0_nlp ...
,dl0_nlp ...
,dl1_nlp ...
,dw0w0_nlp ...
,dw0l0_nlp ...
,dw0l1_nlp ...
,dl0l0_nlp ...
,dl1l1_nlp ...
] = ...
PAD_newton_nlp_XYC_0( ...
 parameter ...
,n_x ...
,n_t ...
,XY_xt__ ...
,omega ...
,l0 ...
,l1 ...
);

str_thisfunction = 'PAD_newton_nlp_XYC_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
nf=0;
n_x = 2; n_t = 128;
XY_xt__ = randn(n_x,n_t);
omega = 0; l0 = 0; l1 = 0;
f_nlp = @(w0l0l1_) PAD_nlp_XYC_strip_0(n_x,n_t,XY_xt__,w0l0l1_(1+0),w0l0l1_(1+1),w0l0l1_(1+2));
[w0l0l1_opt_,fval_opt] = fminsearch(f_nlp,[omega,l0,l1],optimset('MaxFunEvals',1024));
omega = w0l0l1_opt_(1+0); l0 = w0l0l1_opt_(1+1); l1 = w0l0l1_opt_(1+2);
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = 1e-6;
parameter.flag_verbose = 1;
[ ...
 parameter ...
,C_n_xx__ ...
,omega ...
,l0 ...
,l1 ...
,niteration_newton ...
,nlp ...
,dw0_nlp ...
,dl0_nlp ...
,dl1_nlp ...
,dw0w0_nlp ...
,dw0l0_nlp ...
,dw0l1_nlp ...
,dl0l0_nlp ...
,dl1l1_nlp ...
] = ...
PAD_newton_nlp_XYC_0( ...
 parameter ...
,n_x ...
,n_t ...
,XY_xt__ ...
,omega ...
,l0 ...
,l1 ...
);
%%%%%%%%;
n_domega = 48+1;
domega_ = linspace(-pi/2,+pi/2,n_domega);
n_dl1 = 64;
dl1_ = linspace(-2,+2,n_dl1);
[domega_wl__,dl1_wl__] = ndgrid(domega_,dl1_);
nlp_wl__ = zeros(n_domega,n_dl1);
tmp_parameter = parameter; tmp_parameter.flag_verbose = 0;
for ndomega=0:n_domega-1;
for ndl1=0:n_dl1-1;
domega = domega_(1+ndomega);
dl1 = dl1_(1+ndl1);
tmp_omega = omega + domega;
tmp_l0 = l0;
tmp_l1 = l1 + max(l0,l1)*dl1;
[ ...
 tmp_parameter ...
,tmp_nlp ...
] = ...
PAD_nlp_XYC_0( ...
 tmp_parameter ...
,n_x ...
,n_t ...
,XY_xt__ ...
,tmp_omega ...
,tmp_l0 ...
,tmp_l1 ...
);
nlp_wl__(1+ndomega,1+ndl1) = tmp_nlp;
end;%for ndl1=0:n_dl1-1;
end;%for ndomega=0:n_domega-1;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
subplot(1,2,1);
nlp_lim_ = [prctile(nlp_wl__,  0,'all'),prctile(nlp_wl__,100,'all')];
imagesc(nlp_wl__,nlp_lim_); axis image; axisnotick; colorbar;
xlabel('dl1'); ylabel('domega');
subplot(1,2,2);
contourf(dl1_,domega_,nlp_wl__,32); axisnotick;
xlabel('dl1'); ylabel('domega');
%%%%%%%%;
f_nlp = @(w0l0l1_) PAD_nlp_XYC_strip_0(n_x,n_t,XY_xt__,w0l0l1_(1+0),w0l0l1_(1+1),w0l0l1_(1+2));
[w0l0l1_opt_,fval_opt] = fminsearch(f_nlp,[omega,l0,l1],optimset('MaxFunEvals',1024));
omega = w0l0l1_opt_(1+0); l0 = w0l0l1_opt_(1+1); l1 = w0l0l1_opt_(1+2);
%%%%%%%%;
n_domega = 48+1;
domega_ = linspace(-pi/2,+pi/2,n_domega);
n_dl1 = 64;
dl1_ = linspace(-2,+2,n_dl1);
[domega_wl__,dl1_wl__] = ndgrid(domega_,dl1_);
nlp_wl__ = zeros(n_domega,n_dl1);
tmp_parameter = parameter; tmp_parameter.flag_verbose = 0;
for ndomega=0:n_domega-1;
for ndl1=0:n_dl1-1;
domega = domega_(1+ndomega);
dl1 = dl1_(1+ndl1);
tmp_omega = omega + domega;
tmp_l0 = l0;
tmp_l1 = l1 + max(l0,l1)*dl1;
[ ...
 tmp_parameter ...
,tmp_nlp ...
] = ...
PAD_nlp_XYC_0( ...
 tmp_parameter ...
,n_x ...
,n_t ...
,XY_xt__ ...
,tmp_omega ...
,tmp_l0 ...
,tmp_l1 ...
);
nlp_wl__(1+ndomega,1+ndl1) = tmp_nlp;
end;%for ndl1=0:n_dl1-1;
end;%for ndomega=0:n_domega-1;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
subplot(1,2,1);
nlp_lim_ = [prctile(nlp_wl__,  0,'all'),min(max(nlp_lim_),prctile(nlp_wl__,100,'all'))];
imagesc(nlp_wl__,nlp_lim_); axis image; axisnotick; colorbar;
xlabel('dl1'); ylabel('domega');
subplot(1,2,2);
contourf(dl1_,domega_,nlp_wl__,32); axisnotick;
xlabel('dl1'); ylabel('domega');
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); n_t=[]; end; na=na+1;
if (nargin<1+na); XY_xt__=[]; end; na=na+1;
if (nargin<1+na); omega=[]; end; na=na+1;
if (nargin<1+na); l0=[]; end; na=na+1;
if (nargin<1+na); l1=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
if ~isfield(parameter,'n_iteration_newton'); parameter.n_iteration_newton = 128; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
tolerance_master = parameter.tolerance_master;
n_iteration_newton = parameter.n_iteration_newton;
flag_verbose = parameter.flag_verbose;
if isempty(omega); omega = 0; end;
if isempty(l0); l0 = 0; end;
if isempty(l1); l1 = 0; end;

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

niteration_newton = 0;
flag_continue = 1;
tmp_parameter = parameter; tmp_parameter.flag_verbose = 0;
while flag_continue;
omega = mod(omega,2*pi);
[ ...
 tmp_parameter ...
,nlp ...
,dw0_nlp ...
,dl0_nlp ...
,dl1_nlp ...
,dw0w0_nlp ...
,dw0l0_nlp ...
,dw0l1_nlp ...
,dl0l0_nlp ...
,dl1l1_nlp ...
] = ...
PAD_nlp_XYC_0( ...
 tmp_parameter ...
,n_x ...
,n_t ...
,XY_xt__ ...
,omega ...
,l0 ...
,l1 ...
);
dl0w0_nlp = dw0l0_nlp; dl1w0_nlp = dw0l1_nlp; dl0l1_nlp = 0; dl1l0_nlp = 0;
y_3_ = [ omega ; l0 ; l1 ];
f_3_ = [dw0_nlp ; dl0_nlp ; dl1_nlp];
df_33__ = [ ...
	  dw0w0_nlp , dw0l0_nlp , dw0l1_nlp ; ...
	  dl0w0_nlp , dl0l0_nlp , dl0l1_nlp ; ...
	  dl1w0_nlp , dl1l0_nlp , dl1l1_nlp ; ...
          ];
df_pinv_33__ = pinv(df_33__,tolerance_master);
dy_3_ = - df_pinv_33__ * f_3_ ;
if flag_verbose;
disp(sprintf(' %% niteration_newton %d/%d: fnorm(f_3_): %0.16f , fnorm(dy_3_)/max(1e-12,fnorm(y_3_)): %0.16f',niteration_newton,n_iteration_newton,fnorm(f_3_),fnorm(dy_3_)/max(1e-12,fnorm(y_3_))));
end;%if flag_verbose;
flag_continue = (niteration_newton< n_iteration_newton) & (fnorm(dy_3_)/max(1e-12,fnorm(y_3_))>=tolerance_master);
niteration_newton = niteration_newton + 1;
y_3_ = y_3_ + dy_3_;
omega = y_3_(1+0);
l0 = y_3_(1+1);
l1 = y_3_(1+2);
end;%while flag_continue;

[ ...
 tmp_parameter ...
,CtCn__ ...
] = ...
PAD_BtBn_0( ...
 tmp_parameter ...
,omega ...
,l0 ...
,l1 ...
);

C_n_xx__ = sqrtm(CtCn__);
if (flag_verbose); disp(sprintf(' %% C_n_xx__ error: %0.16f',fnorm(CtCn__ - transpose(C_n_xx__)*C_n_xx__)/fnorm(CtCn__))); end;

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;







