function ...
[ ...
 nlp ...
] = ...
PAD_nlp_XYC_strip_0( ...
 n_x ...
,n_t ...
,XY_xt__ ...
,omega ...
,l0 ...
,l1 ...
);

na=0;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); n_t=[]; end; na=na+1;
if (nargin<1+na); XY_xt__=[]; end; na=na+1;
if (nargin<1+na); omega=[]; end; na=na+1;
if (nargin<1+na); l0=[]; end; na=na+1;
if (nargin<1+na); l1=[]; end; na=na+1;

[ ...
 ~ ...
,CtCn__ ...
] = ...
PAD_BtBn_0( ...
 [] ...
,omega ...
,l0 ...
,l1 ...
);

nlp = 0.5 * sum( (CtCn__*XY_xt__).*XY_xt__ , 'all' ) - n_t*( (l0 + l1)/2 - 0.5*n_x*log(2*pi) );







