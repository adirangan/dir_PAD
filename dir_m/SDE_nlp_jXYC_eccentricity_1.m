function ...
[ ...
 nlp ...
] = ...
SDE_nlp_jXYC_eccentricity_1( ...
 n_x ...
,n_j ...
,l0 ...
,l1 ...
);

na=0;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); n_j=[]; end; na=na+1;
if (nargin<1+na); l0=[]; end; na=na+1;
if (nargin<1+na); l1=[]; end; na=na+1;

nlp = 0;
nlp = nlp + n_j*0.5*(l0 - l1).^2; %<-- eccentricity penalty. ;
nlp = nlp + n_j*0.5*(l0.^2 + l1.^2); %<-- magnitude penalty. ;


