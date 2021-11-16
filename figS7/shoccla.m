if (exist('graphics_toolkit') )
 set(0,'defaultaxeslinewidth',1.5)
end
%
f=figure(1, 'position', [100 300 1200 300]);
clf ;
styles={'o-', 's-', 'v-', '*-', '^-'} ; 
colors={'k' 'r' 'g' 'b'} ;
ms=5;
% some parameters
bcmemfr=0.25;

agc1=1;
agc1s= [ 1 1.5 2 ] ;
iagc=0;
for agc1 = agc1s
iagc=iagc+1;

k12=10;
k22=10;
qtcell=0;

nbcrs=[2:5 7 10 15];
occls= [0 0.5 0.9 1] ;

ioc=0;
for occl = occls ;
ioc=ioc+1; % occlusion counter

subplot(1,length(occls),ioc); hold on ; box on ; 

i=0; % counter
leg={};
for dkp = [ 1 10 100 ]
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
 mbc(nbcr)=d.mbcr(end); % ratio of MBCR1/(MBCR1-N) vs simulation # ; 'end' means last simulation
end
leg=[leg {['K_{eq}^{11} / K_{eq}^{i >1,1}=10^{', num2str(0.1*round(10*log10(1000/dkp))),'}']}] ;

plot( nbcrs, mbc(nbcrs), [char(colors(iagc)),char(styles(i))] , 'markersize', ms)
ylim([0 1.])

set(gca, 'fontsize', 14)

end
if (ioc==1)
 l=legend(leg) ; legend boxoff ; 
 set(l, 'fontsize',9)
 ylabel('\it\zeta')
 if (k12==10)
  text(-0.25,1.075,'A','fontsize',17)
 elseif (k12==0)
  text(-0.25,1.075,'B','fontsize',17)
 end
 text(9.5,0.2+iagc*0.1,['\it\alpha_1^T=',num2str(agc1)] ,'fontsize',13, 'color',char(colors(iagc)))
end
title(['Occlusion=', num2str(occl)]);
xlabel('\it #Epitopes');

end

end % ag1s

set(gcf, 'paperpositionmode','auto');
print(gcf, '-dtiff', ['occl-k12=',num2str(k12),'.tif'])
print(gcf, '-depsc2', ['occl-k12=',num2str(k12),'.eps'])
