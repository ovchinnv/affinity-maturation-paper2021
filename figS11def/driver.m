addpath('../') ; format long ;
param ;% basic init
qdrive=1 ;% manual set flag for some vars

occl=1. ; %

nbcr=3;
k12=0 ; agc1=1 ;
k22=10 ; agc2=1. ;
k32=1 ; agc3=1. ;

%%%%%%%%%%%%%% lognormal
ene_=log(affinity);  % generally, the energy values above will have negative values, which are disallowed in the log ; we offset them (arbtrarily)
mu_=ene_(7) ;
sig_=lambda * 0.6 ; % x 0.75 - 1 is a reasonable range

ibc0 = 1/(sig_*sqrt(2*pi)) * exp (-0.5 * ( ( ene_ - mu_ )/sig_).^2 ); % initial affinity dist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize :
ibc0 = ibc0 / sum(ibc0) * minimum_xi * 100; 
ibc1=ibc0 ;
ibc1=ibc0 ;
ibc1=ibc0 ;

run
show