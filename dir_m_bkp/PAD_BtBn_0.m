function ...
[ ...
 parameter ...
,BtBn__ ...
,dw0_BtBn__ ...
,dl0_BtBn__ ...
,dl1_BtBn__ ...
,dw0w0_BtBn__ ...
,dw0l0_BtBn__ ...
,dw0l1_BtBn__ ...
] = ...
PAD_BtBn_0( ...
 parameter ...
,omega ...
,l0 ...
,l1 ...
);

str_thisfunction = 'PAD_BtBn_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
omega = pi/7 ; l0 = randn(); l1 = randn();
domega = pi*1e-5; dl0 = 1e-5; dl1 = 1e-5;
[ ~,BtBn_000__,dw0_BtBn_000__,dl0_BtBn_000__,dl1_BtBn_000__,dw0w0_BtBn_000__,dw0l0_BtBn_000__,dw0l1_BtBn_000__] = PAD_BtBn_0([],omega + 0*domega,l0 + 0*dl0,l1 + 0*dl1);
%%%%;
[ ~,BtBn_p00__,dw0_BtBn_p00__,dl0_BtBn_p00__,dl1_BtBn_p00__,dw0w0_BtBn_p00__,dw0l0_BtBn_p00__,dw0l1_BtBn_p00__] = PAD_BtBn_0([],omega + 1*domega,l0 + 0*dl0,l1 + 0*dl1);
[ ~,BtBn_n00__,dw0_BtBn_n00__,dl0_BtBn_n00__,dl1_BtBn_n00__,dw0w0_BtBn_n00__,dw0l0_BtBn_n00__,dw0l1_BtBn_n00__] = PAD_BtBn_0([],omega - 1*domega,l0 + 0*dl0,l1 + 0*dl1);
%%%%;
[ ~,BtBn_0p0__,dw0_BtBn_0p0__,dl0_BtBn_0p0__,dl1_BtBn_0p0__,dw0w0_BtBn_0p0__,dw0l0_BtBn_0p0__,dw0l1_BtBn_0p0__] = PAD_BtBn_0([],omega + 0*domega,l0 + 1*dl0,l1 + 0*dl1);
[ ~,BtBn_0n0__,dw0_BtBn_0n0__,dl0_BtBn_0n0__,dl1_BtBn_0n0__,dw0w0_BtBn_0n0__,dw0l0_BtBn_0n0__,dw0l1_BtBn_0n0__] = PAD_BtBn_0([],omega + 0*domega,l0 - 1*dl0,l1 + 0*dl1);
%%%%;
[ ~,BtBn_00p__,dw0_BtBn_00p__,dl0_BtBn_00p__,dl1_BtBn_00p__,dw0w0_BtBn_00p__,dw0l0_BtBn_00p__,dw0l1_BtBn_00p__] = PAD_BtBn_0([],omega + 0*domega,l0 + 0*dl0,l1 + 1*dl1);
[ ~,BtBn_00n__,dw0_BtBn_00n__,dl0_BtBn_00n__,dl1_BtBn_00n__,dw0w0_BtBn_00n__,dw0l0_BtBn_00n__,dw0l1_BtBn_00n__] = PAD_BtBn_0([],omega + 0*domega,l0 + 0*dl0,l1 - 1*dl1);
%%%%;
dw0_BtBn_est__ = (BtBn_p00__ - BtBn_n00__)/(2*domega);
dl0_BtBn_est__ = (BtBn_0p0__ - BtBn_0n0__)/(2*dl0);
dl1_BtBn_est__ = (BtBn_00p__ - BtBn_00n__)/(2*dl1);
disp(sprintf(' %% dw0_BtBn_est__ vs dw0_BtBn_000__: %0.16f',fnorm(dw0_BtBn_est__ - dw0_BtBn_000__)/fnorm(dw0_BtBn_est__)));
disp(sprintf(' %% dl0_BtBn_est__ vs dl0_BtBn_000__: %0.16f',fnorm(dl0_BtBn_est__ - dl0_BtBn_000__)/fnorm(dl0_BtBn_est__)));
disp(sprintf(' %% dl1_BtBn_est__ vs dl1_BtBn_000__: %0.16f',fnorm(dl1_BtBn_est__ - dl1_BtBn_000__)/fnorm(dl1_BtBn_est__)));
%%%%;
dw0w0_BtBn_est__ = (dw0_BtBn_p00__ - dw0_BtBn_n00__)/(2*domega);
dw0l0_BtBn_est__ = (dw0_BtBn_0p0__ - dw0_BtBn_0n0__)/(2*dl0);
dw0l1_BtBn_est__ = (dw0_BtBn_00p__ - dw0_BtBn_00n__)/(2*dl1);
disp(sprintf(' %% dw0w0_BtBn_est__ vs dw0w0_BtBn_000__: %0.16f',fnorm(dw0w0_BtBn_est__ - dw0w0_BtBn_000__)/fnorm(dw0w0_BtBn_est__)));
disp(sprintf(' %% dw0l0_BtBn_est__ vs dw0l0_BtBn_000__: %0.16f',fnorm(dw0l0_BtBn_est__ - dw0l0_BtBn_000__)/fnorm(dw0l0_BtBn_est__)));
disp(sprintf(' %% dw0l1_BtBn_est__ vs dw0l1_BtBn_000__: %0.16f',fnorm(dw0l1_BtBn_est__ - dw0l1_BtBn_000__)/fnorm(dw0l1_BtBn_est__)));
%%%%;
dl0l0_BtBn_est__ = (dl0_BtBn_0p0__ - dl0_BtBn_0n0__)/(2*dl0);
disp(sprintf(' %% dl0l0_BtBn_est__ vs dl0l0_BtBn_000__: %0.16f',fnorm(dl0l0_BtBn_est__ -   dl0_BtBn_000__)/fnorm(dl0l0_BtBn_est__)));
dl1l1_BtBn_est__ = (dl1_BtBn_00p__ - dl1_BtBn_00n__)/(2*dl1);
disp(sprintf(' %% dl1l1_BtBn_est__ vs dl1l1_BtBn_000__: %0.16f',fnorm(dl1l1_BtBn_est__ -   dl1_BtBn_000__)/fnorm(dl1l1_BtBn_est__)));
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); omega=[]; end; na=na+1;
if (nargin<1+na); l0=[]; end; na=na+1;
if (nargin<1+na); l1=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
% U_n__ = [ +c -s ; +s +c ];
% BtBn__ = U_n_ * [ expl0 0 ; 0 expl1 ] * transpose(U_n_) ;
%%%%%%%%;
c = cos(omega) ; s = sin(omega) ;
s2w = c*s; c2w = c^2 - s^2; %<-- double angle formulae. ;
expl0 = exp(l0); expl1 = exp(l1);
BtBn__ = [   expl0*c^2 + expl1*s^2   ,   s2w*(expl0-expl1)        ; ...
             s2w*(expl0-expl1)       ,   expl0*s^2 + expl1*c^2   ];

%%%%;
if nargout>2;
%%%%;
% ;
dw0_BtBn__ = [   2*s2w*(expl1 - expl0)  ,    c2w*(expl0 - expl1)      ; ...
		   c2w*(expl0 - expl1)  ,  2*s2w*(expl0 - expl1)     ];
% actually derivative with respect to log(expl0) ;
dl0_BtBn__ = expl0 * [   c^2   ,  +s2w        ; ...
                           +s2w   ,   s^2       ];
% actually derivative with respect to log(expl1) ;
dl1_BtBn__ = expl1 * [   s^2   ,  -s2w        ; ...
                           -s2w   ,   c^2       ];
% ;
dw0w0_BtBn__ = [   2*c2w*(expl1 - expl0)   ,  -4*s2w*(expl0 - expl1)             ; ...
		  -4*s2w*(expl0 - expl1)   ,   2*c2w*(expl0 - expl1)            ];
% ;
dw0l0_BtBn__ = expl0 * [   -2*s2w    ,    +c2w      ; ...
  		                +c2w    ,  +2*s2w     ];
% ;
dw0l1_BtBn__ = expl1 * [   +2*s2w    ,    -c2w      ; ...
  		                -c2w    ,  -2*s2w     ];
%%%%;
end;%if nargout>2;
%%%%;

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;







