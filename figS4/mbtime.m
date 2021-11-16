% plot MBC fraction vs infection number
if (exist('graphics_toolkit') )
% graphics_toolkit('gnuplot');
end
%
f=figure(1, 'position', [100 300 1200 300]);
clf ;
styles={'o', 's', 'v', '*', '^','.'} ; 
% some "constants"
bcmemfr=0.25;
k12=10;
k22=10;
qtcell=0;

nbcrs=[10];
%occls=1 ;
occls = [0 0.5 0.9 1];

ioc=0;
for occl = occls ;
ioc=ioc+1; % occlusion counter

subplot(1,length(occls),ioc); hold on ; box on ; 

leg={};
for nbcr=nbcrs
for iter = 1:50
 fname = [ 'nbc=', num2str(nbcr),...
          '|mfr=', num2str(bcmemfr),...
          '|o=', num2str(occl),...
          '|k12=', num2str(k12),...
          '|k22=', num2str(k22),...
          '|it=', num2str(iter),...
          '|qt=', num2str(qtcell) ] 

 d=load([fname,'.mat']);
 runs=[1:length(d.mbcr)];
 mbc(iter,runs)=d.mbcr(:);
 plot( runs, mbc(iter,runs), ['k',char(styles(end))], 'markersize', 10 )
end


end
% average & std
errorbar( runs, mean(mbc(:,runs),1), std(mbc(:,runs),1), 'ro-' )
%ylim([0 1])
ylim([0 0.6])
set(gca, 'fontsize', 14)

if (ioc==1)
% l=legend(leg) ; legend boxoff ; set(l, 'fontsize',8)
ylabel('\it MBC_{stem}/MBC_{total}');
end
title(['Occlusion=', num2str(occl)]);

xlabel('\it Exposure #');

end

set(gcf, 'paperpositionmode','auto');
print(gcf, '-dtiff', ['mbctime-k12=',num2str(k12),'-nbc=',num2str(nbcrs(1)),'rnd.tif'])

