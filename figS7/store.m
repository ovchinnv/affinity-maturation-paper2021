% save data
% fname is based on simulation parameters
%
fname = [ 'nbc=', num2str(nbcr),...
          '|mfr=', num2str(bcmemfr),...
          '|o=', num2str(occl),...
          '|ag1=', num2str(agc1),...
          '|k12=', num2str(bcr(1).k2),...
          '|k22=', num2str(bcr(2).k2),...
          '|dk1=', num2str(dkpenalty(1)),...
          '|qt=', num2str(qtcell) ] ;
%
save('-mat',[fname,'.mat'], 'mbcr', 'nrun', 'bcmemfr', 'occl', 'nbcrs', 'dkpenalty', 'mbcsave', 'affinity', 'mbcouttot', 'agsout' )

