function ...
[ ...
 nlp ...
] = ...
PAD_nlp_tXYC_eccentricity_0( ...
 n_x ...
,n_t ...
,l0 ...
,l1 ...
);

na=0;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); n_t=[]; end; na=na+1;
if (nargin<1+na); l0=[]; end; na=na+1;
if (nargin<1+na); l1=[]; end; na=na+1;

nlp = 0;
nlp = nlp + n_t*0.5*(l0 - l1).^2; %<-- eccentricity penalty. ;
nlp = nlp + n_t*0.5*(l0.^2 + l1.^2); %<-- magnitude penalty. ;


