% plot MBC fraction vs infection number
if exist('OCTAVE_VERSION')
% graphics_toolkit('gnuplot');
 set(0,'defaultaxeslinewidth',1.5)
end
%
f1=figure(1, 'position', [100 300 1200 300]); clf ;
f2=figure(2, 'position', [100 300 1200 300]); clf ;
styles={'o', 's', 'v', '*', '^','.'} ; 
colors={'k' 'r' 'g' 'b'} ;
ms=5;
% parameters
bcmemfr=0.25;
k12=10;
%k12=0;
k22=10;
qtcell=0;

ags1s={ '1-1-1-1-1-1' } ;

occls = [0 0.5 0.9 1];
dkps=[ 1 10 100];

iag=0;
for ags1 = ags1s
iag=iag+1;

ioc=0;
for occl = occls ;
ioc=ioc+1; % occlusion counter

leg={};
nbcs=[2:5 7 10 15];
nbcs=10;

for nbcr=nbcs
for idkp = 1:length(dkps)
 dkp=dkps(idkp);
 fname = [...,
          'a1-', char(ags1),...
          '|n=', num2str(nbcr),...
          '|mfr=', num2str(bcmemfr),...
          '|o=', num2str(occl),...
          '|k12=', num2str(k12),...
          '|k22=', num2str(k22),...
          '|dk1=', num2str(dkp),...
          '|qt=', num2str(qtcell) ] ;

 d=load([fname,'.mat']);
 runs=[1:length(d.mbcr)];
 mbc(runs)=d.mbcr(:)
 figure(f1);
 subplot(1,length(occls),ioc); hold on ; box on ;
 plot( runs, mbc(runs), [char(colors(iag)), char(styles(idkp)),'-'], 'markersize', ms )
%
 leg=[leg {['K_{eq}^{11} / K_{eq}^{i >1,1}=10^{', num2str(0.1*round(10*log10(1000/dkp))),'}']}] ;

% compute average AB affinity after each run :
 for irun=runs
  mbcs=cell2mat(d.mbcsave(irun)); % 2D array, class vs BCR#
  tas=d.affinity*mbcs ./ sum(mbcs,1);  % average affinity of each AB in this run
  af1(irun)=tas(1);
  af2(irun)=mean(tas(2:end)); % lump the rest into one var
 end % irun
 figure(f2);
 subplot(1,length(occls),ioc); hold on ; box on ;
 plot( runs, af1(runs)./af2(runs), [char(colors(iag)), char(styles(idkp)),'--'], 'markersize', ms )
end

end
% average & std
figure(f1);
subplot(1,length(occls),ioc); hold on ; box on ; 
set(gca, 'fontsize', 14)
set(gca, 'yscale','linear');
ylim([0 1])
if (ioc==1)
 l=legend(leg) ; legend boxoff ; 
 set(l, 'fontsize',12)
 ylabel('\it\zeta');
 if (k12==10)
  text(-0.25,1.075,'A','fontsize',17)
 elseif (k12==0)
  text(-0.25,1.075,'B','fontsize',17)
 endif
 text(1.1,0.2+iag*0.1,['\it\alpha_1^T=',strrep(char(ags1),'-',', ') ] ,'fontsize',13, 'color',char(colors(iag)))

end
title(['Occlusion=', num2str(occl)]);
xlabel('\it Exposure #');
%
figure(f2);
subplot(1,length(occls),ioc); hold on ; box on ; 
set(gca, 'fontsize', 14)
set(gca, 'yscale','log');
ylim([0 3000])

if (ioc==1)

 if (k12==10)
  text(0.,3e3,'C','fontsize',17)
 elseif (k12==0)
  text(0.,3e3,'D','fontsize',17)
 end

 l=legend(leg) ; legend boxoff ; set(l, 'fontsize',12)
 ylabel('\it< \kappa^1 > / < \kappa^{i>1} >', 'interpreter', 'tex');
end
title(['Occlusion=', num2str(occl)]);
xlabel('\it Exposure #');

end

end

figure(f1);
set(gcf, 'paperpositionmode','auto');
print(gcf, '-dtiff', ['mbctime-k12=',num2str(k12),'si.tif']);
print(gcf, '-depsc2', ['mbctime-k12=',num2str(k12),'si.eps']);
figure(f2);
set(gcf, 'paperpositionmode','auto');
print(gcf, '-dtiff', ['afftime-k12=',num2str(k12),'si.tif']);
print(gcf, '-depsc2', ['afftime-k12=',num2str(k12),'si.eps']);
