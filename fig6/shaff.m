if exist('OCTAVE_VERSION')
% graphics_toolkit('gnuplot');
 set(0,'defaultlinelinewidth',1.5)
end
%
f=figure(1, 'position', [100 300 1200 300]);
clf ;
styles={'o-', 's-', 'v-', '*-', '^-'} ; 
% some parameters
bcmemfr=0.25;
agc1=1;
k12=10;
k22=10;
qtcell=0;

nbcrs=[2:5 7 10 15];
occls=[0 0.5 0.9 1] ;

ioc=0;
for occl = occls ;
ioc=ioc+1; % occlusion counter

subplot(1,length(occls),ioc); hold on ; box on ; 

i=0; % counter
leg={};
for dkp = [ 1 10 100 ] ;
i=i+1 ; %counter
for nbcr=nbcrs
 fname = [ 'nbc=', num2str(nbcr),...
          '|mfr=', num2str(bcmemfr),...
          '|o=', num2str(occl),...
          '|ag1=', num2str(agc1),...
          '|k12=', num2str(k12),...
          '|k22=', num2str(k22),...
          '|dk1=', num2str(dkp),...
          '|qt=', num2str(qtcell) ] 

 d=load([fname,'.mat']);
 mbcs=cell2mat(d.mbcsave(end)) ;% MBC output at last run
 tas=d.affinity*mbcs ./ sum(mbcs,1);  % average affinity of each AB in this run
 tar(nbcr) = tas(1) / mean(tas(2:end)) ; % affinity ratio

end
leg=[leg {['K_{11} / K_{1 i > 1}=10^{', num2str(0.1*round(10*log10(1000/dkp))),'}']}] ;

plot( nbcrs, tar(nbcrs), ['k',char(styles(i))] )
set(gca, 'yscale','log')
ylim([0.01 2000])

set(gca, 'fontsize', 14)

end
if (ioc==1)
 l=legend(leg) ; legend boxoff ; 
 set(l, 'fontsize',8)
end
title(['Occlusion=', num2str(occl)]);

xlabel('\it #Epitopes');

end

set(gcf, 'paperpositionmode','auto');
print(gcf, '-dtiff', ['affr-k12=',num2str(k12),'agc1=',num2str(agc1),'.tif'])
print(gcf, '-depsc2', ['affr-k12=',num2str(k12),'agc1=',num2str(agc1),'.eps'])
