function output = shyge(k_,N,K,n_);
%%%%%%%%;
if (nargin<1);
N = 25; K = 15; k_ = -1:K; n_ = k_+3;
for nl=1:length(k_);
ln_1_(nl) = log(hygepdf(k_(nl),N,K,n_(nl)));
end;%for nl=1:length(k_);
ln_2_ = shyge(k_,N,K,n_);
subplot(1,2,1);
plot(k_,ln_1_,'k.-',k_,ln_2_,'ro-');
xlabel('k');ylabel('log(hyge(k,N,K,n))');
subplot(1,2,2);
plot(k_,log10(abs(ln_1_-ln_2_)./abs(ln_1_)),'kx-');
xlabel('k');ylabel('log10(relative error in log(hyge(k,N,K,n))');
disp('returning');return;
end;%if (nargin<1);
%%%%%%%%;
M = N-K; m_ = n_-k_;
output = snchoosek(M,m_) + snchoosek(K,k_) - snchoosek(N,n_);
tmp_ij_ = find(k_<0 | k_>K | m_<0 | m_>M);
output(tmp_ij_) = -Inf;
