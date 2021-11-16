% T-cell model
%
if (qtcell)
% an issue is that is you set tcells too low, the GC dies
 ppag=0.1; % peptides per antigen (note that we really want something like peptides per epitope, since Tcells are epitope-specific)
% empirically, ppag needs to be quite small for Tcells not to get numerous
% it is not easy to justify biologically
%
 Cmhc=1; % this is shorthand for Cmhc = [bcr] * ppag / [mhc] ; e.g. # peptides per mhc
 Tfit = 200*minimum_xi ; % maximum allowed tcr concentration for fit
 mhc2fit=Tfit ; % makes no significant difference
% set minimal and maximal activation (at Tfit):
 htmin=0.01 ; htmax=0.9 ;
 kt0=htmin/((Tfit-mhc2fit*htmin)*(1-htmin)); % binding constant to unloaded MHC at very low B cell number
 kt1=htmax/((Tfit-mhc2fit*htmax)*(1-htmax)); % binding constant to loaded MHC at very low  B cell number
 assert(kt0>=0);
 assert(kt1>=0);
 tcr_tot = 10 * minimum_xi * ones(1,estep); % initial concentration
% note that I use "tcr" but they are actually cells, not receptors
% note that after this point, Tfit does not appear anywhere ; ratio kt1/kt0 independent of Tfit
 tcr_unb = eps * ones(1,estep);
 tcr_min = 0.001 * Tfit ; % minimum number of tcr receptors ; note that this could compromise the GC death mechanism, if too high
% for tfh-cell expansion :
 ktscale = 1 ;
 ktprof = ktscale * 1.25 ; % comparable to parameters in Mayer et al 2019 T cell model
 ktdeath = ktscale * 0.25 ; % if you set the death rate too high, T cells will not expand quickly enough
end
