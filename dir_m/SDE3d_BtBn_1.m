function ...
[ ...
 parameter ...
,BtBn__ ...
] = ...
SDE3d_BtBn_1( ...
 parameter ...
,omega_ ...
,l0 ...
,l1 ...
,l2 ...
);

str_thisfunction = 'SDE3d_BtBn_1';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
rng(0);
omega = pi/7 ; l0 = randn(); l1 = randn();
[ ~,BtBn_0__] = SDE_BtBn_0([],omega,l0,l1);
[ ~,BtBn_1__] = SDE3d_BtBn_1([],omega,l0,l1);
disp(sprintf(' %% 2d: BtBn_0__ vs BtBn_1__: %0.16f',fnorm(BtBn_0__-BtBn_1__)/fnorm(BtBn_0__)));
omega_ = 2*pi*rand(3,1); l2 = randn();
[ ~,BtBn_2__] = SDE3d_BtBn_1([],omega_,l0,l1,l2);
L0L1L2_ = exp([l0;l1;l2]); L0L1L2_ = sort(L0L1L2_,'descend');
[tmp_] = svds(BtBn_2__,3); tmp_ = sort(tmp_,'descend');
disp(sprintf(' %% 3d: L0L1L2_ vs tmp_: %0.16f',fnorm(L0L1L2_ - tmp_)/fnorm(L0L1L2_)));
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); omega_=[]; end; na=na+1;
if (nargin<1+na); l0=[]; end; na=na+1;
if (nargin<1+na); l1=[]; end; na=na+1;
if (nargin<1+na); l2=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(omega_); omega_=0; end;
if isempty(l0); l0=0; end;
if isempty(l1); l1=0; end;
if isempty(l2); l2=0; end;

if ( (numel(omega_)==1) & (numel(l0)==1) & (numel(l1)==1) ); 
%%%%%%%%;
% interpret omega as angle of dominant principal-component in 2-dimensions. ;
% and (l0,l1) as the log-principal-values. ;
%%%%%%%%;
[parameter,BtBn__] = SDE_BtBn_0(parameter,omega_,l0,l1);
else;
%%%%%%%%;
% interpret omega_ as the euler-angles for an orthonormal basis in 3-dimensions. ;
% and (l0,l1,l2) as the log-principal-values. ;
%%%%%%%%;
assert(numel(omega_)==3);
eazimub = omega_(1+0);
epolara = omega_(1+1);
egammaz = omega_(1+2);
%%%%%%%%;
% The general formula used here (and in get_template_1.m) is as follows. ;
% let sa and ca be sin(polar_a) and cos(polar_a), respectively. ;
% let sb and cb be sin(azimu_b) and cos(azimu_b), respectively. ;
% let sc and cc be sin(gamma_z) and cos(gamma_z), respectively. ;
% And rotation by azimu_b about the +z-axis is represented as: ;
% Rz(azimu_b) = ;
% [ +cb -sb 0 ] ;
% [ +sb +cb 0 ] ;
% [  0   0  1 ] ;
% And rotation by polar_a about the +y-axis is represented as: ;
% Ry(polar_a) = ;
% [ +ca 0 +sa ] ;
% [  0  1  0  ] ;
% [ -sa 0 +ca ] ;
% And rotation by gamma_z about the +z-axis is represented as: ;
% Rz(gamma_z) = ;
% [ +cc -sc 0 ] ;
% [ +sc +cc 0 ] ;
% [  0   0  1 ] ;
% Which, collectively, implies that under the transform: ;
% Rz(azimu_b) * Ry(polar_a) * Rz(gamma_z), ;
% Which is the same as: ;
% [ +cb -sb 0 ] [ +ca*cc -ca*sc +sa ]   [ +cb*ca*cc - sb*sc , -cb*ca*sc -sb*cc , +cb*sa ];
% [ +sb +cb 0 ] [ +sc    +cc    0   ] = [ +sb*ca*cc + cb*sc , -sb*ca*sc +cb*cc , +sb*sa ];
% [  0   0  1 ] [ -sa*cc +sa*sc +ca ]   [ -sa*cc            , +sa*sc           , +ca    ];
%%%%%%%%;
cb = cos(eazimub); sb = sin(eazimub);
ca = cos(epolara); sa = sin(epolara);
cc = cos(egammaz); sc = sin(egammaz);
U_n__= [ +cb*ca*cc - sb*sc , -cb*ca*sc - sb*cc , +cb*sa ; ...
         +sb*ca*cc + cb*sc , -sb*ca*sc + cb*cc , +sb*sa ; ...
         -sa*cc            , +sa*sc           , +ca    ];
expl0 = exp(l0); expl1 = exp(l1); expl2 = exp(l2);
BtBn__ = U_n__*diag([expl0;expl1;expl2])*transpose(U_n__);
%%%%%%%%;
end;%if ( (numel(omega_)==1) & (numel(l0)==1) & (numel(l1)==1) ); 

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;







