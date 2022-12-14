function [r_pv__,r_rtn__,r_mu__,r_sg__,r_sum__,r_cap__] = cumulative_intersect_ver4(r_remaining_A_,rdrop_A_,r_remaining_B_,rdrop_B_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% function cumulative_intersect_ver4(r_remaining_A_,rdrop_A_,r_remaining_B_,rdrop_B_);
% ;
% This function calculates the cumulative intersection of rdrop_A_ and rdrop_B_. ;
% This function *strongly* assumes that rdrop_A_ and rdrop_B_ contain exactly the same set of indices. ;
% Both these index-sets must range from 1 to n_total_A==n_total_B. ;
% Each of these index-sets must also be unique. ;
% ;
% Inputs: ;
% r_remaining_A_       integer array of size n_iteration_A-by-1 holding the number of rows remaining at the start of each iteration (i.e., r_remaining_A_(1)-r_remaining_A_(2) = number of rows eliminated during first iteration). ;
% rdrop_A_           integer array of size n_total_A-by-1 listing the row indices in the order they are eliminated (i.e., rdrop_A_(end) is retained the longest). ;
% r_remaining_B_       integer array of size n_iteration_B-by-1 holding the number of rows remaining at the start of each iteration (i.e., r_remaining_B_(1)-r_remaining_B_(2) = number of rows eliminated during first iteration). ;
% rdrop_B_           integer array of size n_total_B-by-1 listing the row indices in the order they are eliminated (i.e., rdrop_B_(end) is retained the longest). ;
%
% Outputs: ;
% r_pv__   double array of size n_iteration_A-by-n_iteration_B. r_pv__(j,k) is the p-value associated with each entry in r_rtn__. ;
% r_rtn__  integer array of size n_iteration_A-by-n_iteration_B. r_rtn__(j,k) is the size of the intersection between rdrop_A and rdrop_B at the beginning of iteration j (for trial-A) and iteration k (for trial-B). ;
% r_mu__   double array of size n_iteration_A-by-n_iteration_B. r_mu__(j,k) is the expected value associated with r_rtn__ assuming random removal. ;
% r_sg__   double array of size n_iteration_A-by-n_iteration_B. r_mu__(j,k) is the standard-deviation of r_rtn__ assuming random removal. ;
% r_sum__  integer array of size n_iteration_A-by-n_iteration_B. r_sum__(j,k) is the (double) cumulative sum of r_cap__. i.e., the total number of ejected rows which lie in the intersection of steps 1-through-j (trial A) and steps 1-through-k (trial B). ;
% r_cap__  integer array of size n_iteration_A-by-n_iteration_B. r_cap__(j,k) is the length of the intersection of the ejected rows from step j (trial A) and step k (trial B). ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (nargin<4);
disp(sprintf(' %% testing cumulative_intersect_ver4'));
rng(1);
rdrop_A_ = transpose(randperm(9));%rdrop_A_ = [9;1;8;2;7;3;6;4;5]; 
r_remaining_A_ = [9;7;5;3;3;2;1;0;0]; 
rdrop_B_ = transpose(randperm(9));%rdrop_B_ = [2;4;6;8;1;3;5;7;9];
r_remaining_B_ = [9;6;6;4;2];
[r_pv__,r_rtn__,r_mu__,r_sg__,r_sum__,r_cap__] = cumulative_intersect_ver4(r_remaining_A_,rdrop_A_,r_remaining_B_,rdrop_B_);
%disp(sprintf(' %% r_cap__: ')); disp(full(r_cap__));
%disp(sprintf(' %% r_sum__: ')); disp(full(r_sum__));
disp(sprintf(' %% r_rtn__: ')); disp(full(r_rtn__));
%disp(sprintf(' %% r_pv__: ')); disp(full(r_pv__));
%disp(sprintf(' %% r_mu__: ')); disp(full(r_mu__));
%disp(sprintf(' %% r_sg__: ')); disp(full(r_sg__));
r_chk__ = zeros(length(r_remaining_A_),length(r_remaining_B_));
for nrA=1:length(r_remaining_A_); for nrB=1:length(r_remaining_B_);
r_chk__(nrA,nrB) = length(intersect(rdrop_A_(end-r_remaining_A_(nrA)+1:end),rdrop_B_(end-r_remaining_B_(nrB)+1:end)));
end;end;%for nrA=1:length(r_remaining_A_); for nrB=1:length(r_remaining_B_);
disp(sprintf(' %% r_chk__: ')); disp(full(r_chk__));
disp(sprintf(' %% error %d ',sum(sum(abs(r_rtn__-r_chk__)))));
disp('returning');return;
end;%if (nargin<4);

%%%%%%%%;
n_iteration_A = length(r_remaining_A_); %<-- total number of iterations. ;
n_total_A = length(rdrop_A_); assert(n_total_A==r_remaining_A_(1)); %<-- total number of rows. ;
r_removed_A_ = -diff([r_remaining_A_;0]);
%%%%%%%%;
n_iteration_B = length(r_remaining_B_); %<-- total number of iterations. ;
n_total_B = length(rdrop_B_); assert(n_total_B==r_remaining_B_(1)); %<-- total number of rows. ;
r_removed_B_ = -diff([r_remaining_B_;0]);
%%%%%%%%;
if (n_total_B~=n_total_A); disp(sprintf(' %% Warning! n_total_A %d ~= n_total_B %d in cumulative_intersect_ver4',n_total_A,n_total_B)); end;
n_total = n_total_A;
%%%%%%%%;
[~,iA_] = sort(rdrop_A_); 
[~,iB_] = sort(rdrop_B_); 
iC__ = sparse(iA_,iB_,1,n_total_A,n_total_B);
%%%%%%%%;
tmp_r_A = zeros(n_total_A,1);
nr_A=0;
for ni_A=1:n_iteration_A;
if (nr_A<n_total_A); tmp_r_A(1+nr_A) = tmp_r_A(1+nr_A) + 1; end;
nr_A = nr_A + r_removed_A_(ni_A);
end;%for ni_A=1:n_iteration_A;
tmp_r_A = cumsum(tmp_r_A);
tmp_r_A_ = sparse(1:n_total_A,tmp_r_A,1,n_total_A,n_iteration_A);
%%%%%%%%;
tmp_r_B = zeros(n_total_B,1);
nr_B=0;
for ni_B=1:n_iteration_B;
if (nr_B<n_total_B); tmp_r_B(1+nr_B) = tmp_r_B(1+nr_B) + 1; end;
nr_B = nr_B + r_removed_B_(ni_B);
end;%for ni_B=1:n_iteration_B;
tmp_r_B = cumsum(tmp_r_B);
tmp_r_B_ = sparse(1:n_total_B,tmp_r_B,1,n_total_B,n_iteration_B);
%%%%%%%%;
r_cap__ = transpose(tmp_r_A_)*iC__*tmp_r_B_; %<-- r_cap__(j,k) is the intersection of the ejected rows from step j (trial A) and step k (trial B). ;
r_sum__ = cumsum(cumsum(r_cap__,1),2); %<-- r_sum__(j,k) is the (double) cumulative sum of r_cap__. i.e., the total number of ejected rows which lie in the intersection of steps 1-j (trial A) and steps 1-k (trial B). ;
r_rmv_A__ = repmat(cumsum([0;r_removed_A_]),1,1+n_iteration_B);
r_rmv_B__ = repmat(transpose(cumsum([0;r_removed_B_])),1+n_iteration_A,1);
r_rtn__ = n_total - r_rmv_A__ - r_rmv_B__;
r_rtn__(2:end,2:end) = r_rtn__(2:end,2:end) + r_sum__; 
r_rtn__ = r_rtn__(1:end-1,1:end-1); %<-- r_rtn__(j,k) is the number of rows remaining at the start of step j (trial A) and step k (trial B). ;
%%%%%%%%;
r_pA__ = (n_total_A - r_rmv_A__)/n_total_A;
r_pB__ = (n_total_B - r_rmv_B__)/n_total_B;
r_mu__ = (r_pA__.*r_pB__) * n_total;
r_sg__ = (r_pA__.*r_pB__) .* (1 - r_pA__.*r_pB__) * n_total;
r_mu__ = r_mu__(1:end-1,1:end-1); r_sg__ = r_sg__(1:end-1,1:end-1);
r_pv__ = 0.5*( 1 - erf((r_rtn__ - r_mu__)./(sqrt(2) * r_sg__)) );
r_pv__(find(~isfinite(r_pv__)))=0.5;
