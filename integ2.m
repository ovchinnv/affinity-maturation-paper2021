% integrate by one step using explicit Euler
%
nriter=0; % keep track of NR iterations
a=zeros(nag,1);
f=zeros(nag,1);
for iag=1:nag
 a(iag,1)=ags(iag).alpha(istep);
end
%
do
fp=zeros(nag);
for iag=1:nag % precompute bound fractions & derivatives
 f(iag,1)=a(iag)-ags(iag).alphatot(istep) ; % first part of constraint
 fp(iag,iag)=1.0 ;
 hx(iag)=0;
 hpx(iag)=0;
 for ibcr=ags(iag).bcr; % allowing multiple antibodies to bind this antigen (though each AB can still only bind one!)
  p=bcr(ibcr).affinity;%*a(iag);
  q=p*bcr(ibcr).k2*a(iag) ;
  d=1./(1+a(iag)*(p+q)) ;
  h=(p+2*q).*d ; % need x a below
  hp=h.*(1-a(iag)*h) + 2*q.*d ; % derivative (h') for NR iteration
  h=a(iag)*h;
% save some values
  bcr(ibcr).hc(istep,:)=h;
  xcorr=bcr(ibcr).xc(istep,:);
  if (qab) % to include AB competition in conservation equation
   xcorr=xcorr + bcr(ibcr).xab(istep,:) ;
   if (qmkmem) % ABs derived from memory cells
    xcorr=xcorr + Cmemab * bcr(ibcr).xmab(istep,:) ;
   end
  end % accounts for competition with Abs
  hx(iag) = hx(iag) + h*xcorr';
  hpx(iag) = hpx(iag) + hp*xcorr';
 end
end
%
for iag=1:nag
 for jag=1:nag
  f(iag) = f(iag) + hx(jag) * occlusion(iag,jag);
  fp(iag,jag) = fp(iag,jag) + hpx(jag) * occlusion(iag,jag); % derivative (f') for NR
 end %jag
end
a = a - fp\f; % NR iteration
%return
nriter=nriter+1;
until (norm(f)<1e-6)
% update arrays :
for iag=1:nag
 ags(iag).alpha(istep:istep+1)=a(iag); % also assign guess for next iteration
% also update total antigen :
 ags(iag).alphatot(istep+1) = ( 1.0 - dt * ags(iag).kalpha_death ) * ags(iag).alphatot(istep) ; % explicit euler
% semi-implicit C-N : more accurate, but it does not matter here because the cell & AB dynamics are dominated by BCR & AB
%alphatot(istep) = ( 2.0 - dt * kalpha_death ) / ( 2.0 + dt * kalpha_death ) * alphatot(istep-1) ; 
end

%%%%%%%%%%%%%%%%%%%%%%%% alpha & therefore hc are known, proceed to integrate %%%%%%%
% account for maximum GC size via modified kprof :
xall=[bcr.xc](istep,:); % concatenate xc's side-by-side
gcsize=sum(xall);
k2_prof=k_prof*(1.0 - gcsize/gcmax);
%
% T-cell model :
if (qtcell)
 h_mhc = min(1, Cmhc*([bcr.hc](istep,:))/(1+qavid)); % MHC activation factor or fraction MHC bound
% i.e. Cmhc converts h_bcr to h_mhc, subject to a maximum activation cutoff
 mhct = ppag / Cmhc * gcsize ; %total MHC (note that this is the only role for ppag in the model)
% mhcb = ppag / Cmhc * sum(xall .* h_mhc ) ;% peptide-bound MHC
% equiv :
 mhcb = mhct * sum(xall .* h_mhc ) / gcsize ;% can roughly view h_mhc as bcr fraction converted to mhc fraction
% solve nondimensional version via kp93
 xt=([mhct-mhcb mhcb]/tcr_tot(istep))'; % unloaded & loaded receptors, resp.
 talpha=tcr_unb(istep)/tcr_tot(istep); %initial guess for NR
 taffinity=[kt0 kt1]*tcr_tot(istep); % make nondimensional per PK93
 do
  ht = 1 - 1./(1.+talpha*taffinity);
  tf = talpha + ht * xt - 1. ; % TCR conservation condition
  tfpal =  talpha + (ht.*(1.- ht))*xt; % f prime times alpha
  talpha = talpha - tf*talpha/tfpal; % NR iteration to compute unbound tcr
%  nriter=nriter+1;
 until abs(tf)<1e-4
 hactivation = ( ht(1) + h_mhc * ( ht(2) - ht(1) ) );
 tcr_unb(istep:istep+1)=tcr_tot(istep)*talpha ; % include initial guess for next iteration
% save hactivation for each bcr
 hactivation=reshape(hactivation,nclass,[]);
 for ibcr=1:nbcr
  bcr(ibcr).hct(istep,:) = hactivation(:,ibcr);
 end
%===== t cell evolution =====%
 tcr_tot(istep+1) = tcr_tot(istep) + dt * ( ktprof * ( tcr_tot(istep)-tcr_unb(istep) ) - ktdeath * max(0,tcr_tot(istep)-tcr_min) ) ;
end
%
for ibcr=1:nbcr % over all bcells

xcorr = bcr(ibcr).xc(istep,:);
%
% adjust proliferation based on maximum allowed clone size (if set)
if (exist('clonemax'))
 k3_prof = k2_prof * (1.0 - sum(xcorr)/clonemax);
else
 k3_prof = k2_prof ;
end
%
% compute dxdt in steps:
dxdt = - k_death * xcorr ; % apoptosis
%
if (qminxi)                        % whether to zero out very small concentrations (not used in the paper)
 xcorr=xcorr.*(xcorr>minimum_xi) ; % note that this would cause dxdt=0, but does not by itself "zero out" x ;
end
%
if (qtcell)
 hactivation=bcr(ibcr).hct(istep,:) ;
else
 hactivation=bcr(ibcr).hc(istep,:)/(1+qavid) ; % binding fraction as activation ; simplest (KP93) model
end
% apply scaling & exponent (see methods section in paper)
hactivation = hcscale * ( hactivation.^hcexp ) ;
%
hxcorr = hactivation.*xcorr ;
dxdt = dxdt + k_death * hxcorr; % deactivated apoptosis (i.e survival signal)

if (~qhprof) % turn off activated proliferation
 hxcorr=xcorr ;
end
%
if (rate_arrays) % rate array form
 dxdt = dxdt + k3_prof.*hxcorr ;% taking + k_prof rather than -k_prof means we can omit 1 in m_ii ; diff from mtrans below
% split the remaining terms into interior and boundary points :
% interior
 dxdt(2:end-1) = dxdt(2:end-1) + 2*k3_prof * ( hxcorr(1:end-2).*bcr(ibcr).w(2:end-1) + ... 
                                               hxcorr(2:end-1).*bcr(ibcr).o(2:end-1) + ... 
                                               hxcorr(3:end).*  bcr(ibcr).e(2:bcr(ibcr).nclass-1) ) ;
% boundary
 dxdt(1)   = dxdt(1)   + 2 * k3_prof * ( hxcorr(1)  *bcr(ibcr).o(1)   + hxcorr(2)*bcr(ibcr).e(1) );
 dxdt(end) = dxdt(end) + 2 * k3_prof * ( hxcorr(end)*bcr(ibcr).o(end) + hxcorr(end-1)*bcr(ibcr).w(end) );
%
else % matrix form
 dxdt = dxdt - hxcorr * k3_prof  +  2*k3_prof * hxcorr * bcr(ibcr).mtrans ; % Eq. 9
% NOTE : cannot factor hxcorr above because constants k_ will get promoted to matrices and added to mtrans and give wrong answer
% a mistake I made in the past
end % rate_arrays
% save dxdt
bcr(ibcr).dxc(istep,:) = dxdt ;
bcr(ibcr).xc(istep+1,:) = bcr(ibcr).xc(istep,:) + dt * dxdt ; % Euler step
%
%== Extensions : plasma and memory cells
%
if (qplasma || qmkmem)
% plasma cell activation function (could be 1 for p=q=0 ; param.m & paper methods)
% note that I am using the same activation, which is not strictly necessary
 hex=pqnorm * (hactivation.^pexit) .* ((1 - hactivation).^qexit); % cell exiting proliferation/mutation cycle
 dexdt = (1-hactivation) .* hex .* xcorr ;
%
 if (qplasma)
  bcr(ibcr).xpl(istep+1,:) = (1 - dt*kpl_death) * bcr(ibcr).xpl(istep,:) + dt * Cpl * dexdt ;
  if (qab)
   bcr(ibcr).xab(istep+1,:) = (1 - dt*kab_death) * bcr(ibcr).xab(istep,:) + dt * kab_prod * bcr(ibcr).xpl(istep,:);
  endif
 end
%
 if (qmkmem)
  bcr(ibcr).xmem(istep+1,:) = (1 - dt*kmem_death) * bcr(ibcr).xmem(istep,:) + dt * Cmem * dexdt ;
  if (qplasma) % 5/1 another source of plasma cells -- see Methods eq. 14
   bcr(ibcr).xmpl(istep+1,:) = (1 - dt*kpl_death) * bcr(ibcr).xmpl(istep,:) + dt * kmpl_prod * bcr(ibcr).xmem(istep,:) ;
   if (qab) % AB secretion by PL derived from memory
    bcr(ibcr).xmab(istep+1,:) = (1 - dt*kab_death) * bcr(ibcr).xmab(istep,:) + dt * kab_prod * bcr(ibcr).xmpl(istep,:);
   endif
  end
 end
%
end

end % loop over ibcr
istep=istep + 1; % increment step counter
