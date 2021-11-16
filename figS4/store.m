% save data
% fname is based on simulation parameters
%
fname = [ 'nbc=', num2str(nbcr),...
          '|mfr=', num2str(bcmemfr),...
          '|o=', num2str(occl),...
          '|k12=', num2str(bcr(1).k2),...
          '|k22=', num2str(bcr(2).k2),...
          '|it=', num2str(iter),...
          '|qt=', num2str(qtcell) ] ;
%
save('-mat',[fname,'.mat'], 'mbcr', 'nrun', 'bcmemfr', 'occl', 'nbcrs', 'dkpenalty', 'dkps', 'mbcouttot', 'agsout' )

