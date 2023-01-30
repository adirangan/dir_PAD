%%%%%%%%;
% Quick start for matlab use. ;
%%%%%%%%;
Notes: Matlab uses 1-based indexing by default, ;
but in this README I'll use zero-based indexing. ;
So an array A_ of size N=4 would have entries: ;
A_(1), A_(2), A_(3), A_(4), ;
which I'll reference by:
A_(1+na), with na ranging from 0:N-1. ;

The main driver is this routine:
SDE_nlp_ijXaABYC_update_all_1.m
with version number indicated at the end (i.e., '1').

The calling sequence for this routine looks like:

function ...
[ ...
 parameter ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_est_ixt___ ...
,n_a ...
,a_est_xa__ ...
,A_est_xx__ ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
] = ...
SDE_nlp_ijXaABYC_update_all_1( ...
 parameter ...
,n_i ...
,n_t_i_ ...
,t_it__ ...
,n_x ...
,X_est_ixt___ ...
,n_a ...
,a_est_xa__ ...
,A_est_xx__ ...
,B_est_omega ...
,B_est_l0 ...
,B_est_l1 ...
,n_j_i_ ...
,index_nt_from_nj_i__ ...
,ignore_Y_tru_ixj___ ...
,Y_tru_ixj___ ...
,C_est_omega ...
,C_est_l0 ...
,C_est_l1 ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Inputs and Outputs: ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

parameter : matlab structure. Holds extra user specified parameters (described below). ;
n_i : integer. number of individuals (e.g., individual animals) recorded. ;
n_t_i_ : integer array of size (n_i,1). number of unique time-points recorded for individual ni. ;
       e.g., individual ni has n_t_i_(1+ni) time-points associated with it. ;
t_it__ : cell array of size (n_i,1). ;
       after defining t_t_ = t_it__{1+ni} and n_t = n_t_i_(1+ni), ;
       t_t_ will be a double array of size (n_t,1). ;
       	    sorted (ascending) time-points associated with individual ni. ;
n_x : integer. number of variables considered (must be 2). ;
X_est_ixt___ : cell array of size (n_i,1). ;
	     after defining t_t_ = t_it__{1+ni} and n_t = n_t_i_(1+ni), ;
             X_est_xt__ will be a double array of size (n_x,n_t). ;
       		  Contains the estimated X-values for individual ni. ;
		  X_est_xt__(1+nx,1+nt) contains the measured variable nx at time nt. ;
n_a : integer. 1+degree of polynomial associated with the drift (i.e., velocity) term in the SDE. ;
    For example, n_a = 1 implies a constant. ;
    n_a = 2 implies a linear velocity term. ;
    n_a = 3 implies a quadratic velocity term. ;
a_est_xa__ : double array of size (n_x,n_a). ;
	   estimated velocity term for the SDE. ;
	   a_est_xa__(1+nx,:) are the polynomial coefficients for variable nx. ;
	   The ordering is from lowest to highest degree. ;
	   So a_est_xa__(1+nx,1+0) is the constant term. ;
	   So a_est_xa__(1+nx,1+1) is the linear term. ;
	   So a_est_xa__(1+nx,1+2) is the quadratic term, and so forth. ;
	   Note that this is the 'fliplr' of the matlab polyval ordering. ;
A_est_xx__ : double array of size (n_x,n_x). ;
	   estimated deterministic interaction matrix for the SDE. ;
%%%%%%%%;
% We assume that the stochastic interaction matrix B ;
% is represented using the eigendecomposition of the inverse-covariance matrix:
% inv(B*transpose(B)) = U_n__ * [exp(l0) 0 ; 0 exp(l1) ] * transpose(U_n__). ;
% with U_n__ = [ +c -s ; +s +c ] ;
% and c = cos(omega), s = sin(omega). ;
%%%%;
B_est_omega : double. omega (above). ;
B_est_l0 : double. l0 (above). ;
B_est_l1 : double. l1 (above). ;
%%%%%%%%;
n_j_i_ : integer array. size (n_i,1). number of measurements for individual ni. ;
index_nt_from_nj_i__ : cell array of size (ni,1). ;
		     after defining t_t_ = t_it__{1+ni}. ;
		     and n_j = n_j_i_(1+ni). ;
		     then index_nt_from_nj_ = index_nt_from_nj_i__{1+ni}. ;
		     will be an integer array of size (n_j,1). ;
		     	  index_nt_from_nj_ can be used to cross-reference measurements to times. ;
			  So, for example, if we look at measurement nj from individual ni, ;
			  the corresponding time (of that measurement) will be: ;
			  t_t_(1+index_nt_from_nj_(1+ni)). ;
ignore_Y_tru_ixj___ : cell array of size (ni,1). ;
		    after defining: ;
		    n_j = n_j_i_(1+ni) and ;
		    ignore_Y_tru_xj__ = ignore_Y_tru_ixj___{1+ni}, ;
		    then ignore_Y_tru_xj__ will be a integer array of size (nx,n_j) ;
		    containing a flag indicating whether or not to ignore the ;
		    measurements for individual ni. ;
		    	       Specifically, ignore_Y_tru_xj__(1+nx,1+nj) ;
			       will be equal to 0 if the measurement at that nx and nj ;
			       is to be used, and will be equal to 1 if the measurement
			       is to be ignored (e.g., due to missing data). ;
Y_tru_ixj___ : cell array of size (ni,1). ;
		    after defining: ;
		    n_j = n_j_i_(1+ni) and ;
		    Y_tru_xj__ = Y_tru_ixj___{1+ni}, ;
		    then Y_tru_xj__ will be a double array of size (nx,n_j) ;
		    containing the actual measurements for individual ni. ;
		    	       Specifically, Y_tru_xj__(1+nx,1+nj) ;
			       will be the variable nx at measurement nj for individual ni.
%%%%%%%%;
% We assume that the stochastic interaction matrix C ;
% is represented using the eigendecomposition of the inverse-covariance matrix:
% inv(C*transpose(C)) = U_n__ * [exp(l0) 0 ; 0 exp(l1) ] * transpose(U_n__). ;
% with U_n__ = [ +c -s ; +s +c ] ;
% and c = cos(omega), s = sin(omega). ;
%%%%;
C_est_omega : double. omega (above). ;
C_est_l0 : double. l0 (above). ;
C_est_l1 : double. l1 (above). ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Usage: ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

Despite the daunting number of inputs listed above, ;
most can actually be ignored. ;
The only inputs that are required are: ;
n_i ;
n_t_i_ ;
t_it__ ;
n_j_i_ ;
index_nt_from_nj_i__ ;
ignore_Y_tru_ixj___ ;
Y_tru_ixj___ ;
which specify how many individuals are measured, ;
the times at which they are measured, ;
the number of measurements for each individual, ;
and which times those measurements correspond to ;
(as well as which measurements should be ignored). ;

We also require: ;
n_x ;
n_a ;
where n_x is the dimension of the system (so, 2). ;
and n_a determines the degree of the polynomial term for the velocity (e.g., 1 for a constant velocity term). ;

The other inputs, such as: ;
X_est_ixt___ ;
a_est_xa__ ;
A_est_xx__ ;
B_est_omega ;
B_est_l0 ;
B_est_l1 ;
C_est_omega ;
C_est_l0 ;
C_est_l1 ;
can all be passed in as empty (i.e., []). ;

The routine will assume reasonable defaults:
X_est_ixt___ zero ;
a_est_xa__ zero ;
A_est_xx__ zero ;
B_est_omega zero ;
B_est_l0 zero ;
B_est_l1 zero ;
C_est_omega zero ;
C_est_l0 zero ;
C_est_l1 zero ;

This means that the default settings for X, a and A are all zero, ;
while B and C are assumed to be the identity matrix. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% parameter settings. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
The parameter input should be a matlab structure, defined like: ;
parameter = struct('type','parameter');
with additional fields defined via: ;
parameter.tolerance_master = 1e-2 ;
parameter.flag_verbose = 1 ;
etcetera.

The fields used are:
parameter.tolerance_master : double. master tolerance used to estimate convergence. default: 1e-2;
parameter.flag_verbose : integer. level of verbosity when running. default: 0;
parameter.flag_disp : integer. flag indicating whether or not to plot. default: 0;
parameter.flag_regularize_eccentricity_BtBn : double. prefactor associated with the regularization term when initially optimizing B and C. default 1;
parameter.flag_regularize_eccentricity_simultaneous : double. prefactor associated with the regularization term during main optimization. default: 0;
parameter.n_iteration_BtBn : integer. maximum number of iterations when initially optimizing B and C. default: 16;
parameter.MaxFunEvals_use_BtBn : integer. maximum number of function-evaluations used in nelder-meade when initially optimizing B and C. default: 1024;
parameter.MaxFunEvals_use_simultaneous : integer. maximum number of function-evaluations used in nelder-meade during main optimization. default: 1024;
