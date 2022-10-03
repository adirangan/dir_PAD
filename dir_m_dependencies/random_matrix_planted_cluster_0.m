function [A_n_,label_A_,n_label_A_,pf_,pi_,snr_,MX_,NX_,mu_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);

verbose=0;

continue_flag=1;
while (continue_flag);
if n_cluster==1; X_ = 0.5; end;
if n_cluster>1; X_ = linspace(0.3,0.7,n_cluster); end;
MX_ = ceil(M.^X_); NX_ = ceil(N.^X_); 
if (sum(MX_)>M | sum(NX_)>N); n_cluster = n_cluster-1; continue_flag=1; end;
if (sum(MX_)<=M & sum(NX_)<=N); continue_flag=0; end;
end;%while (continue_flag);

%%%%%%%%;
% Here the signal-to-noise-ratio (snr) is defined to be the ratio between: ;
% numerator: the dominant singular-value of the signal, and ;
% denominator: the dominant singular-value of the noise. ;
% We can estimate this as follows: ;
% Given a matrix A of size MA-by-NA of the form: ;
%  A_n = [ B_n C_n ] ;
%        [ D_n E_n ] ;
% with each entry drawn from N(0,1), ;
% except for entries of B which are drawn from N(0,1) + mu, ;
% we see that the covariance matrix is equal to (in expectation): ;
% A_n*A_t = [ B_n*B_t + C_n*C_t , B_n*D_t + C_n*E_t ] ;
%           [ D_n*B_t + E_n*C_t , D_n*D_t + E_n*E_t ] ;
%         = [ NB*mu^2*1_n*1_t + NA*I , 0    ] ;
%           [ 0                      , NA*I ] ;
% where 1_n*1_t is the MB-by-MB ones-matrix. ;
% The dominant eigenvector of the 'signal' is [ 1_n ; 0 ]/sqrt(MB), ;
% with dominant eigenvalue  MB*NB*mu^2, ;
% and the dominant eigenvalue of the 'noise' is NA. ;
% Consequently, we define the snr to be: ;
% snr = mu*sqrt(MB*NB) / sqrt(NA). ;
% Note that this is derived using the covariance-matrix A_n*A_t, ;
% rather than A_t*A_n, which is what we expect when calculating ;
% correlations between rows (i.e., searching for row-clusters). ;
%%%%%%%%;
mu_ = sqrt(N)*snr./sqrt(MX_.*NX_);
snr_ = sqrt(MX_.*NX_).*mu_/sqrt(N);

A_n_ = randn(M,N); 
for nj=1:2; 
pf_{nj} = randperm(size(A_n_,nj)); 
[~,pi_{nj}] = sort(pf_{nj}); [~,pi_{nj}] = sort(pi_{nj}); 
end;% for nj=1:2;
B_n_ = cell(n_cluster,1);
B_label_ = zeros(M,1); 
MX_sum = 0; NX_sum = 0;
for ncluster=1:n_cluster;
MX = MX_(ncluster); NX = NX_(ncluster); mu = mu_(ncluster);
if (verbose>0); disp(sprintf(' %% ncluster %d: snr %0.2f; X %0.2f; M %d MX %d; N %d NX %d mu %0.2f',ncluster,snr,X_(ncluster),M,MX,N,NX,mu)); end;
e_n_ = (2*(rand(NX,1)>0.5) - 1)/sqrt(NX); 
B_n_{ncluster}=randn(MX,NX)+mu*ones(MX,1)*transpose(e_n_)*sqrt(NX);
if (verbose>0); disp(sprintf(' %% B_n_*e_n_: mean %0.2f std %0.2f',mean(B_n_{ncluster}*e_n_),std(B_n_{ncluster}*e_n_))); end;
MX_ij_ = MX_sum + (1:MX); MX_c_ij_ = setdiff(1:M,MX_ij_);
NX_ij_ = NX_sum + (1:NX); NX_c_ij_ = setdiff(1:N,NX_ij_);
A_n_(pf_{1}(MX_ij_),pf_{2}(NX_ij_)) = B_n_{ncluster};
if (verbose>0); 
tmp_A_n_ = A_n_(pi_{1}(MX_c_ij_),pi_{2}(NX_ij_));
disp(sprintf(' %% A_n_(pi_{1}(MX_c_ij_),pi_{2}(NX_ij_))*e_n_: mean %0.2f std %0.2f',mean(tmp_A_n_*e_n_),std(tmp_A_n_*e_n_)));
end;%if (verbose>0); 
B_label_(pi_{1}(MX_sum + (1:MX))) = ncluster;
MX_sum = MX_sum + MX; NX_sum = NX_sum + NX;
end;%for ncluster=1:n_cluster;
[label_A_,n_label_A_] = label_num_to_enum_0(B_label_);