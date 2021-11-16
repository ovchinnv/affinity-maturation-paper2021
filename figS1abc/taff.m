%
% =================================================================
% plot Perelson's total (or average) affinity measure (A)
ifig=ifig+1;figure(ifig);clf;hold on;box on;%grid on;
leg={};
for iag=1:nag
 for ibcr=ags(iag).bcr
  ta=(bcr(ibcr).xc*bcr(ibcr).affinity')./sum(bcr(ibcr).xc,2) ; % average affinity instead :
%  ta=(bcr(ibcr).xc*bcr(ibcr).affinity')./minimum_xi ;
  plot(time,ta,[char(colors(ibcr))] ) ;
  keq=eval(['k',num2str(ibcr),'2']);
  leg=[leg { ['BCR #',num2str(ibcr),', K^2_{eq} = ',num2str(keq),' \times K^1_{eq}'] } ];
%  leg=[leg { ['AG #',num2str(iag),'; BCR #',num2str(ibcr)] } ];
 end
end
% affinity of other cells :
%
if (1*qmkmem) % memory cells
 for iag=1:nag
  for ibcr=ags(iag).bcr
   ta=(bcr(ibcr).xmem*bcr(ibcr).affinity')./sum(bcr(ibcr).xmem,2) ;
   plot(time,ta,[char(colors(ibcr)),'-.'] ) ;
   keq=eval(['k',num2str(ibcr),'2']);
   leg=[leg { ['MBC #',num2str(ibcr),', K^2_{eq} = ',num2str(keq),' \times K^1_{eq}'] } ];
%   leg=[leg { ['AG #',num2str(iag),'; MBC #',num2str(ibcr)] } ];
  end
 end
end

set(gca, 'yscale','log', 'fontsize',16);
%l=legend(leg, 'location', 'southeast'); 
l=legend(leg, 'location', 'northwest'); 
%set(l, 'orientation', 'horizontal');
legend boxoff ; 
set(l, 'fontsize', 11)
ylabel('\it Average affinity (\Sigma \kappa_i x_i/\Sigma x_i)', 'fontsize', 14) ;
xlabel('\it t(days)', 'fontsize', 14) ;
xlim([0 max_t]);
%
ylim([0 1e5]);
if (occl<0.6)
 text(-5.5,1.5e5,'C','fontsize',17)
elseif (occl>=0.9)
 text(-5.5,1.5e5,'F','fontsize',17)
end
set(gcf, 'paperpositionmode','auto')
print(gcf, '-depsc2', 'A.eps')
print(gcf, '-dpng', 'A.png')
%
