if exist('OCTAVE_VERSION')
 set(0,'defaultlinelinewidth',1.5)
end

ifig=0; % figure index

qexp=1 ; % whether to show exp

markers={'s','o','*','x','v','^','d'};
%colors={'g', 'r', 'b', 'm', 'c', 'k', 'g', 'r', 'b', 'm'};

time = dt * [0:nstep];

% plot total xc
ifig=ifig+1;figure(ifig) ; clf ; hold on ; box on ;
leg={};
for iag=1:nag
 for ibcr=ags(iag).bcr
  tx=sum(bcr(ibcr).xc,2)/minimum_xi ;
%  plot(time,tx,[char(colors(ibcr))] ) ;
  keq=eval(['k',num2str(ibcr),'2']);
  leg=[leg { ['BCR #',num2str(ibcr),', K^2_{eq} = ',num2str(keq),' \times K^1_{eq}'] } ];
 end
end
%
if (qexp)
 gcexp ;
 if (nag==1); leg={'Simulation'};end
 leg=[leg {'Wittenbrink et al. 2011'}];
end
l=legend(leg, 'location', 'northeast'); legend boxoff ; set(l, 'fontsize', 12)
set(gca, 'yscale','linear', 'fontsize',14);
ylim([0 5e3]);
xlim([0 35]);
text(-6,5200,'A','fontsize',17)
%text(-6,5200,'D','fontsize',17)
ylabel('\it GC clone size (cells)', 'fontsize', 14) ;
xlabel('\it t(days)', 'fontsize', 14) ;
print(gcf, '-dpng', 'gcsize')
print(gcf, '-depsc2', 'gcsize')

gcmbc
if (qexp)
 if (nag==1); leg={'Simulation'};end
 leg=[leg {'Weisel et al. 2016'}];
end

l=legend(leg, 'location', 'northeast'); 
set(gca, 'yscale','linear', 'fontsize',14);
legend boxoff ; set(l, 'fontsize', 12)
text(-4,3.75,'B','fontsize',17)
%=============================================
% output MBC for inset
axmain=gca;
if (1*qmkmem)
 leg=[];
% ifig=ifig+1; figure(ifig); clf ; hold on ; box on ; %, 'position', [ 100 100 675 420 ]) ; %clf
 for iag=1:nag
  for ibcr=ags(iag).bcr
   tmem=sum(bcr(ibcr).xmem,2)/minimum_xi ;
   tmems(ibcr)=tmem(end);
   leg=[leg { ['MBC#',num2str(ibcr)] } ];
  end
 end

 if (nag>1)
% b=bar(iag,tmems(iag)) ;
 aginset=axes('position', [0.64,0.45,0.25,0.25]); hold on ; box on;
 for iag=1:nag
  b=bar(iag,tmems(iag),0.33)
  set(b,'facecolor',char(colors(iag))) ; 
%  set(b, 'facecolor', ['r', 'g', 'b'])
 end
 set(gca,'xtick',[1:length(leg)])
 set(gca,'xticklabel',leg)
 set(gca,'fontsize',9, 'fontangle','italic', 'fontweight','bold')
 end
end
set(gcf,'currentaxes', axmain)
%=============================================

print(gcf, '-dpng', 'dmbc')
print(gcf, '-depsc2', 'dmbc')

if (qexp)
 gcplc
 if (nag==1); leg={'Simulation'};end
 leg=[leg {'Weisel et al. 2016'}];
end
%
l=legend(leg, 'location', 'southeast'); 
%l=legend(leg, 'location', 'northeast'); 
legend boxoff ; set(l, 'fontsize', 12)
text(-5,10.5,'C','fontsize',17)
set(gca, 'yscale','linear','fontsize', 14);
print(gcf, '-dpng', 'dplc')
print(gcf, '-depsc2', 'dplc')

