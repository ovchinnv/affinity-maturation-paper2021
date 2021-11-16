if (exist('OCTAVE_VERSION'))
% graphics_toolkit('gnuplot');
 graphics_toolkit('fltk');
% graphics_toolkit('qt');
end
% draw fit ; optimization of h

load 'fit-07-Nov-2021.mat' ;

figure(1) ;  hold on ; box on ;
plot(rr,'k-*');
plot(xpar(:,1),'r-o');
plot(xpar(:,2),'g-s');

xlabel('\it iteration')

set(gca, 'fontsize', 17)

legend('\it M.S.Err.', '\it\epsilon', '\it\rho','location', 'northeastoutside'); 
legend boxoff ; 

set(gcf, 'paperpositionmode','auto')
print(gcf, '-depsc2', 'hfit.eps')
print(gcf, '-dpng', 'hfit.png')
