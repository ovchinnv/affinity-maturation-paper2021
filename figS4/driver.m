addpath('../') ; format long ;
% to do basic init :
param
qdrive=1 ;% manual set flag for some vars
%
% perform several runs in one script
%
if (~exist('occlusions'))
 occlusions=[ 0.9 0.99 0.999 0.9999 1] ;% see is passed in
end
%
%
if (~exist('nbcs'))
 nbcs=[ 2:5 7 10 15 ] ;% see is passed in
end

for occl = occlusions ;
for iter = iters
for nbc=nbcs

nrun=6 ; % number of vaccinations/infections
nbcrs = nbc * ones(1, nrun) ; % BCR/epitope pairs for each run

if (~exist('bcmemfr'))
 bcmemfr = 0.25 ; % fraction of GC seed that is taken from memory cells
end
%
%
if (~exist('k12'))
 k12=10 ;
end

ibc0=zeros(1,nclass);ibc0(8)=minimum_xi ; % simple initial affinity dist
ibc1=ibc0 ;

clear mbcr mbcouttot agsout;
dkps=10.^( sqrt(0.5)*randn(1,nrun) + 1 )  % sample from gaussian in log space mean advantage over other Abs : 10^1 = 10

for irun=1:nrun
 nbcr=nbcrs(irun)
% remaining initial conditions (now that nbcr is set for this run):
 for ibcr=2:nbcr
  eval(['k',num2str(ibcr),'2=10;']);
  if (irun==1) % init only for run 1; subsequently, hybrid init from memory
   eval(['ibc',num2str(ibcr),'=ibc0;']);
  end
 end
 run
% show mbc ratio
 mbcout=reshape([bcr.xmem](istep,:), nclass, [] ); % memory concentration at the end
 mbcouttot(irun,:)=sum(mbcout) ; % sum over classes
 mbcr(irun)=mbcouttot(irun,1)/sum(mbcouttot(irun,1:end)) % ratio of desired (stem) to other mbcs
% antigens
 agsout(irun,:)=reshape([ags.alpha],[],nag)(istep,:) % save final antigen concentration
% now, reinitialize for next run ; take a fraction of naive cells ; the remaining cells are from the memory pool
 dkpenalty=1000*ones(1,nbcr) ; % default memory cell affinity penalty on reinfection
 dkp=dkps(irun) ;
 dkpenalty(1)=dkp; 
% but disadvantaged according to assumed ag drift ; normalize to the same initial conc, and rerun GC reaction 
 for ibcr=1:nbcr % go over all bcrs and compute new initial conditions
  aff=bcr(ibcr).affinity ;
  newaff=aff/dkpenalty(ibcr) ; % penalize memory cells by affinity specified for each bcr
  ibc = interp1(log(newaff), mbcout(:,ibcr), log(aff), 'extrap', 0 ) ; % unknown values are 0
% normalize :
  ibc = ibc/sum(ibc) * sum(ibc0) * bcmemfr + (1-bcmemfr) * ibc0 ;
  sum(ibc-ibc0) % make sure this is zero
  eval ( ['ibc', num2str(ibcr),'=ibc;' ] ) ; % set bc in a script
%
 end
end
store
end % nbc
end % dkps or iter
end % occlusions
