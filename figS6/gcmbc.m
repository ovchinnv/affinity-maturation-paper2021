% compute mbc production rate from simulation data
% from integ2 :
% hex=pqnorm * (hactivation.^pexit) .* ((1 - hactivation).^qexit); % cell exiting proliferation/mutation cycle
% dexdt = (1-hactivation) .* hex .* xcorr ;
% bcr(ibcr).xmem(istep+1,:) = (1 - dt*kmem_death) * bcr(ibcr).xmem(istep,:) + dt * Cmem * dexdt ;
if (~exist('qmkmem'))
 qmkmem=0;
end
if (qmkmem==0)
 return
end

ifig=ifig+1;figure(ifig)
leg=[];
clf ; hold on ; box on ;
scale=1/(75) ;
% to use finite difference with existing MBC counts:
for iag=1:nag
 for ibcr=ags(iag).bcr
   tmem=sum(bcr(ibcr).xmem,2);
   difftmem=sum(diff(bcr(ibcr).xmem),2);
   dmemdt=(difftmem/dt + kmem_death * tmem(1:end-1))/minimum_xi ;
   plot(time(1:end-1),scale*dmemdt,[char(colors(ibcr)),'-'] ) ;
   keq=eval(['k',num2str(ibcr),'2']);
   leg=[leg { ['MBC #',num2str(ibcr),', K^2_{eq} = ',num2str(keq),' \times K^1_{eq}'] } ];
 end
end

ylabel('\it MBC production (a.u.)', 'fontsize', 14) ;
xlabel('\it t(days)', 'fontsize', 14) ;
xlim([0 35]);

