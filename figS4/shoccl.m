if exist('OCTAVE_VERSION')
 set(0,'defaultaxeslinewidth',1.5)
 graphics_toolkit('gnuplot');
end
%
f=figure(1, 'position', [100 300 1200 300]);
clf ;
styles={'o', 's', 'v', '*', '^','.'} ; 
% set some parameters
bcmemfr=0.25;
k12=0; % monovalent
%k12=20; % bivalent
k22=10;
qtcell=0;

nbcrs=[2:5 7 10 15];
occls = [ 0 0.5 0.9 1];

ioc=0;
for occl = occls ;
ioc=ioc+1; % occlusion counter

subplot(1,length(occls),ioc); hold on ; box on ; 

i=0; % counter
leg={};
for iter = 1:50
 i=i+1 ; %counter
for nbcr=nbcrs
 fname = [ 'nbc=', num2str(nbcr),...
          '|mfr=', num2str(bcmemfr),...
          '|o=', num2str(occl),...
          '|k12=', num2str(k12),...
          '|k22=', num2str(k22),...
          '|it=', num2str(iter),...
          '|qt=', num2str(qtcell) ] 

 d=load([fname,'.mat']);
 mbc(iter,nbcr)=d.mbcr(end);
end
leg=[leg {['i=',num2str(iter)] } ];

plot( nbcrs, mbc(iter,nbcrs), ['k',char(styles(end))], 'markersize', 10 )
%set(gca, 'yscale','log')
%set(gca, 'xscale','log')
ylim([0 0.8])
xlim([2 16])

set(gca, 'fontsize', 14)

end
% average & std
errorbar( nbcrs, mean(mbc(:,nbcrs),1), std(mbc(:,nbcrs),1), 'r.-' )

if (ioc==1)
% l=legend(leg) ; legend boxoff ; set(l, 'fontsize',8)
% ylabel('\it MBC_{stem}/MBC_{total}');
 ylabel('\it\zeta')
% text(-2,.87,'B','fontsize',17)
 if (k12==10)
  text(-2,.87,'A','fontsize',17)
  axmain=gca;
  axpdf=axes('position', [0.1625, 0.49, 0.1, 0.37])
  dkpen;
  set(gca, 'fontsize',8, 'ytick',[]) ; %box off ; 
  set(l, 'fontsize',14) ;% increase xlabel size
  set(gca, 'linewidth', 1)
  set(gcf, 'currentaxes',axmain);
 else
  text(-2,.87,'B','fontsize',17)
 end
end
title(['Occlusion=', num2str(occl)]);

xlabel('\it #Epitopes');

end

set(gcf, 'paperpositionmode','auto');
print(gcf, '-dtiff', ['occl-k12=',num2str(k12),'rnd.tif'])
print(gcf, '-depsc2', ['occl-k12=',num2str(k12),'rnd.eps'])

