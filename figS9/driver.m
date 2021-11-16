addpath('../') ; format long ;
param ;% basic init
qdrive=1 ;% manual set flag for some vars
occl=0. ; %
nbcr=1;
k12=1 ; agc1=1 ;
%%%%%%%%%%%%%% lognormal
ene_=log(affinity);  % generally, the energy values above will have negative values, which are disallowed in log ; we offset them (arbtrarily)
mu_=ene_(7) ;
sig_=lambda * 0.6 ; % x 0.75 - 1 is a reasonable range
%
ibc0 = 1/(sig_*sqrt(2*pi)) * exp (-0.5 * ( ( ene_ - mu_ )/sig_).^2 );  % initial affinity dist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize :
ibc0 = ibc0 / sum(ibc0) * minimum_xi * 100;
ibc1 = ibc0;
run
