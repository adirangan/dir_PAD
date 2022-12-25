function [pu,lpu] = push_pull_block_pu_0(N,f_pos,f_neg,V_0,V_1,n_pos,n_neg);
%%%%%%%%;
% \begin{eqnarray} ;
% p_{p}\leq p_{u} = \frac{N!}{(N-V-V')!\cdot V!\cdot V'!} \left[ \sum_{n=n_{+}}^{n=VV'} \binom{VV'}{n} \cdot f_{+}^{n}\cdot(1-f_{+})^{VV'-n} \right] \left[ \sum_{n=n_{-}}^{n=VV'} \binom{VV'}{n} \cdot f_{-}^{n}\cdot(1-f_{-})^{VV'-n} \right] ;
% \label{Eq_pu} ;
% \end{eqnarray} ;
%%%%%%%%;

W = V_0 * V_1;

pf_pos = 0;
for n=n_pos:W;
pf_pos = pf_pos + exp( lnchoosek(W,n) + n*log(f_pos) + (W-n)*log(1-f_pos) ) ;
end;%for n=n_pos:W;
lpf_pos = log(pf_pos);

pf_neg = 0;
for n=n_neg:W;
pf_neg = pf_neg + exp( lnchoosek(W,n) + n*log(f_neg) + (W-n)*log(1-f_neg) ) ;
end;%for n=n_neg:W;
lpf_neg = log(pf_neg);

lpu = lnchoosek(N,V_0) + lnchoosek(N-V_0,V_1) + lpf_pos + lpf_neg ;

pu = exp(lpu);

