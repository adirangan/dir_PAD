function [r_pv__,r_rtn__,r_mu__,r_sg__,r_sum__,r_cap__] = cumulative_intersect_ver5(r_remaining_A_,rdrop_A_,r_remaining_B_,rdrop_B_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% function cumulative_intersect_ver5(r_remaining_A_,rdrop_A_,r_remaining_B_,rdrop_B_);
% ;
% This function calculates the cumulative intersection of rdrop_A_ and rdrop_B_. ;
% This function does not assume that rdrop_A_ and rdrop_B_ contain exactly the same set of indices. ;
% Moreover, this function does not assume that these index sets are unique; rdrop_A_ and rdrop_B_ can have repeating indices. ;
% The value of r_rtn__(j,k) is the length of the intersection between trial-A and trial-B at the beginning of iteration j (for trial-A) and k (for trial-B). ;
% This intersection is calculated as follows: sum over all unique indices, adding the minimum of mA and mB for each index, ;
% where mA and mB are the multiplicities of that index in trial-A and trial-B, respectively. ;
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
disp(sprintf(' %% testing cumulative_intersect_ver5'));
for nrng=1:100;
rng(nrng); 
n_A = 18;
rdrop_A_ = round(8*rand(n_A,1)+2);
r_remaining_A_ = [18;14;9;7;5;3;2;1]; 
n_B = 21;
rdrop_B_ = round(8*rand(n_B,1)+4);
r_remaining_B_ = [21;18;16;12;9;6;2];
[r_pv__,r_rtn__,r_mu__,r_sg__,r_sum__,r_cap__] = cumulative_intersect_ver5(r_remaining_A_,rdrop_A_,r_remaining_B_,rdrop_B_);
%disp(sprintf(' %% r_cap__: ')); disp(full(r_cap__));
%disp(sprintf(' %% r_sum__: ')); disp(full(r_sum__));
%disp(sprintf(' %% r_rtn__: ')); disp(full(r_rtn__));
%disp(sprintf(' %% r_pv__: ')); disp(full(r_pv__));
%disp(sprintf(' %% r_mu__: ')); disp(full(r_mu__));
%disp(sprintf(' %% r_sg__: ')); disp(full(r_sg__));
r_chk__ = zeros(length(r_remaining_A_),length(r_remaining_B_));
%%%%%%%%;
for nrA=1:length(r_remaining_A_); for nrB=1:length(r_remaining_B_);
%%%%%%%%;
ao_ = rdrop_A_(1+end-r_remaining_A_(nrA):end);
bo_ = rdrop_B_(1+end-r_remaining_B_(nrB):end);
[au_,ao_to_au_,au_to_ao_] = unique(ao_);
ao_by_au_xref_ = sparse(1:length(ao_),au_to_ao_,1,length(ao_),length(au_));
[bu_,bo_to_bu_,bu_to_bo_] = unique(bo_);
bo_by_bu_xref_ = sparse(1:length(bo_),bu_to_bo_,1,length(bo_),length(bu_));
[cu_,au_to_cu_,bu_to_cu_] = intersect(au_,bu_,'stable'); %<-- cu_ stores unique labels in the intersection of au_ and bu_. ;
au_by_cu_xref_ = sparse(au_to_cu_,1:length(cu_),1,length(au_),length(cu_));
bu_by_cu_xref_ = sparse(bu_to_cu_,1:length(cu_),1,length(bu_),length(cu_));
%%%%%%%%;
m_a_ = zeros(length(cu_),1); %<-- multiplicity of elements of unique intersection within a. ;
m_b_ = zeros(length(cu_),1); %<-- multiplicity of elements of unique intersection within b. ;
for ncu=1:length(cu_);
nau = find(sum(au_by_cu_xref_(:,ncu),2)); assert(length(nau)==1);
tmp_ao_ = find(sum(ao_by_au_xref_(:,nau),2)); assert(length(tmp_ao_)>=1);
m_a_(ncu) = length(tmp_ao_);
nbu = sum(find(bu_by_cu_xref_(:,ncu),2)); assert(length(nbu)==1);
tmp_bo_ = find(sum(bo_by_bu_xref_(:,nbu),2)); assert(length(tmp_bo_)>=1);
m_b_(ncu) = length(tmp_bo_);
end;%for ncu=1:length(cu_);
%%%%%%%%;
m_x_ = min(m_a_,m_b_); %<-- minimum of multiplicities. ;
n_m = sum(m_x_); %<-- size of intersection. ;
r_chk__(nrA,nrB) = n_m;
%%%%%%%%;
end;end;%for nrA=1:length(r_remaining_A_); for nrB=1:length(r_remaining_B_);
%%%%%%%%;
%disp(sprintf(' %% r_chk__: ')); disp(full(r_chk__));
disp(sprintf(' %% nrng %d error %d ',nrng,sum(sum(abs(r_rtn__-r_chk__)))));
end;%for nrng=1:100;
disp('returning');return;
end;%if (nargin<4);

ao_ = rdrop_A_; bo_ = rdrop_B_;
[au_,ao_to_au_,au_to_ao_] = unique(ao_);
ao_by_au_xref_ = sparse(1:length(ao_),au_to_ao_,1,length(ao_),length(au_));
[bu_,bo_to_bu_,bu_to_bo_] = unique(bo_);
bo_by_bu_xref_ = sparse(1:length(bo_),bu_to_bo_,1,length(bo_),length(bu_));
[cu_,au_to_cu_,bu_to_cu_] = intersect(au_,bu_,'stable'); %<-- cu_ stores unique labels in the intersection of au_ and bu_. ;
au_by_cu_xref_ = sparse(au_to_cu_,1:length(cu_),1,length(au_),length(cu_));
bu_by_cu_xref_ = sparse(bu_to_cu_,1:length(cu_),1,length(bu_),length(cu_));

m_a_ = zeros(length(cu_),1); %<-- multiplicity of elements of unique intersection within a. ;
m_b_ = zeros(length(cu_),1); %<-- multiplicity of elements of unique intersection within b. ;
for ncu=1:length(cu_);
nau = find(sum(au_by_cu_xref_(:,ncu),2)); assert(length(nau)==1);
tmp_ao_ = find(sum(ao_by_au_xref_(:,nau),2)); assert(length(tmp_ao_)>=1);
m_a_(ncu) = length(tmp_ao_);
nbu = sum(find(bu_by_cu_xref_(:,ncu),2)); assert(length(nbu)==1);
tmp_bo_ = find(sum(bo_by_bu_xref_(:,nbu),2)); assert(length(tmp_bo_)>=1);
m_b_(ncu) = length(tmp_bo_);
end;%for ncu=1:length(cu_);

m_x_ = min(m_a_,m_b_); %<-- minimum of multiplicities. ;
n_x_ = cumsum([0;m_x_]); %<-- one less than lowest label for unique element cu_(ncu) ;
al_ = zeros(size(ao_)); bl_ = zeros(size(bo_)); %<-- labels. ;
for ncu=1:length(cu_);
tmp_m = m_x_(ncu);
nau = find(sum(au_by_cu_xref_(:,ncu),2)); assert(length(nau)==1);
tmp_ao_ = find(sum(ao_by_au_xref_(:,nau),2)); assert(length(tmp_ao_)>=tmp_m);
tmp_as_ = sort(tmp_ao_,'descend');
al_(tmp_as_(1:tmp_m)) = n_x_(ncu) + (1:tmp_m);
nbu = sum(find(bu_by_cu_xref_(:,ncu),2)); assert(length(nbu)==1);
tmp_bo_ = find(sum(bo_by_bu_xref_(:,nbu),2)); assert(length(tmp_bo_)>=tmp_m);
tmp_bs_ = sort(tmp_bo_,'descend');
bl_(tmp_bs_(1:tmp_m)) = n_x_(ncu) + (1:tmp_m);
end;%for ncu=1:length(cu_);

a0_ = zeros(size(ao_)); a0_(find(al_==0))=1; a0c_ = flip(cumsum(a0_,'reverse')); ra_ = r_remaining_A_ - a0c_(r_remaining_A_); %<-- here we calculate the labels remaining (per-iteration) in trial-A after removing unused labels. ;
b0_ = zeros(size(bo_)); b0_(find(bl_==0))=1; b0c_ = flip(cumsum(b0_,'reverse')); rb_ = r_remaining_B_ - b0c_(r_remaining_B_); %<-- here we calculate the labels remaining (per-iteration) in trial-B after removing unused labels. ;
ax_ = al_(find(al_>0)); %<-- now we remove unused labels. ;
bx_ = bl_(find(bl_>0)); %<-- now we remove unused labels. ;

%size(ra_),; size(ax_),; size(rb_),; size(bx_),;
[r_pv__,r_rtn__,r_mu__,r_sg__,r_sum__,r_cap__] = cumulative_intersect_ver4(ra_,ax_,rb_,bx_);


