function ...
[ ...
 nlp ...
] = ...
SDE_nlp_dtZPB_strip_0( ...
 n_dt ...
,dt_dt_ ...
,n_x ...
,ZP_xdt__ ...
,B_omega ...
,B_l0 ...
,B_l1 ...
);

na=0;
if (nargin<1+na); n_dt=[]; end; na=na+1;
if (nargin<1+na); dt_dt_=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); ZP_xdt__=[]; end; na=na+1;
if (nargin<1+na); B_omega=[]; end; na=na+1;
if (nargin<1+na); B_l0=[]; end; na=na+1;
if (nargin<1+na); B_l1=[]; end; na=na+1;

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

nlp = 0.5 * sum( sum( (BtBn_xx__*ZP_xdt__).*ZP_xdt__ , 1 ) ./ reshape(max(1e-12,dt_dt_),[1,n_dt]) , 2 ) - n_dt*( (B_l0 + B_l1)/2 - 0.5*n_x*log(2*pi) ) + sum( 0.5*n_x*log(dt_dt_) ) ;







