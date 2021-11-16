%
 param ;
%
 if (qmkmem)
  if ~exist('Cmemab')
   Cmemab=0; % antibody production from memory cells ; does not belong here, but included in case missing from param.m, for compat.
  end
  %
  if(qplasma)
   if ~exist('kmpl_prod')
    kmpl_prod=0 ; % rate of PL production from memory cells ; icluded in case missing from param.m
   end
  end
 end
%
 t=it; % set initial time
%
 time0=now;

 while istep<estep
  integ2 ;
  t=t+dt;
  if (mod(t,t_out)<dt)
   fprintf('%s%12.5f\n', ' Time : ',t)
  end
 end
% hc(i) is computed at i+1 iteration because it is needed for x(i+1); so therefore we do not have h(nstep) here;
% for plotting/analysis purposes, populate it :
 for ibcr=1:nbcr
  bcr(ibcr).hc(istep,:)=bcr(ibcr).hc(istep-1,:);
  if (qtcell)
   bcr(ibcr).hct(istep,:)=bcr(ibcr).hct(istep-1,:);
  end
 end
%
 time1=now;
 sec_per_day = 60 * 60 * 24 ;
 dtime=(time1-time0) * sec_per_day
