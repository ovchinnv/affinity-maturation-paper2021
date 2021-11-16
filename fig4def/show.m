if exist('OCTAVE_VERSION')
 set(0,'defaultlinelinewidth',1.5)
end

ifig=0; % figure index

qexp=0 ; % whether to show previous data

markers={'s','o','*','x','v','^','d'};
colors={'r', 'g', 'b', 'm', 'c', 'k', 'r', 'g', 'b', 'm'};
colors={'b', 'g', 'r', 'm', 'c', 'k', 'g', 'r', 'b', 'm'};

time = dt * [0:nstep];

taff ;

% plot total xc
ifig=ifig+1;figure(ifig) ; clf ; hold on ; box on ;
leg={};
for iag=1:nag
 for ibcr=ags(iag).bcr
  tx=sum(bcr(ibcr).xc,2)/minimum_xi ;
  plot(time,tx,[char(colors(ibcr))] ) ;
%  leg=[leg { ['AG #',num2str(iag),'; BCR #',num2str(ibcr)] } ];
  keq=eval(['k',num2str(ibcr),'2']);
  leg=[leg { ['BCR #',num2str(ibcr),', K^2_{eq} = ',num2str(keq),' \times K^1_{eq}'] } ];
 end
end

ylabel('\it GC size (cells)', 'fontsize', 14) ;
xlabel('\it t(days)', 'fontsize', 14) ;
xlim([0 max_t]);
l=legend(leg, 'location', 'northeast'); legend boxoff ; set(l, 'fontsize', 12)
set(gca, 'yscale','linear', 'fontsize',14);
ylim([0 5e3]);
if (occl==0)
 text(-6,5200,'A','fontsize',17)
elseif (occl==1)
 text(-6,5200,'D','fontsize',17)
end
%
print(gcf, '-dpng', 'gcsize')
print(gcf, '-depsc2', 'gcsize')

gcmbc

l=legend(leg, 'location', 'northeast'); 
set(gca, 'yscale','linear', 'fontsize',14);
legend boxoff ; set(l, 'fontsize', 10)
if (occl==0)
 text(-4.5,2.6,'B','fontsize',17)
elseif (occl==1)
 text(-4.5,2.6,'E','fontsize',17)
end
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
% print(gcf, '-dpng', 'membar')
% print(gcf, '-depsc2', 'membar')
 end
end
set(gcf,'currentaxes', axmain)
%=============================================
print(gcf, '-dpng', 'dmbc')
print(gcf, '-depsc2', 'dmbc')

