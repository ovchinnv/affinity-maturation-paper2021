% compute mbc production rate from simulation data
% from integ2 :
% hex=pqnorm * (hactivation.^pexit) .* ((1 - hactivation).^qexit); % cell exiting proliferation/mutation cycle
% dexdt = (1-hactivation) .* hex .* xcorr ;
% bcr(ibcr).xmem(istep+1,:) = (1 - dt*kmem_death) * bcr(ibcr).xmem(istep,:) + dt * Cmem * dexdt ;
if (~exist('ifig')) ; ifig=2 ; end
ifig=ifig+1;figure(ifig)
leg=[];
clf ; hold on ; box on ;
scale=1/(75) ;
% to use finite difference with existing PC counts:
for iag=1:nag
 for ibcr=ags(iag).bcr
   pl=bcr(ibcr).xpl + bcr(ibcr).xmpl;
%   pl=bcr(ibcr).xmpl; % these cells _might_ be relevant ; they are derived from memry cells, and could be of low enough aff;
%    pl=bcr(ibcr).pl; % these cells peak early, then decay, because they are closely derived from GC B cells
   tpl=sum(pl,2);
   difftpl=sum(diff(pl),2);
   dpldt=(difftpl/dt + kpl_death * tpl(1:end-1))/minimum_xi ; % exclude death
%   dmemdt=(difftpl/dt)/minimum_xi ;
   plot(time(1:end-1),scale*dpldt,[char(colors(ibcr)),'-'] ) ;
   leg=[leg { ['AG #',num2str(iag),'; PLC #',num2str(ibcr)] } ];
 end
end
% plot exp data
d=load('../data/pelissier4f.dat');

t=d(1:3:end,1);
gc=d(1:3:end,2);
gcu=d(2:3:end,2)-gc;
gcl=gc-d(3:3:end,2);
% crude scale of minimum & maximum values to get estimate of standard deviation :
% to get the SD estimate (sig), assume that the data are Gaussian distributed and use the estimate :
% x_max ~ 0.9 * sig * sqrt( 2 * log (nsamp) ) ; where nsamp is the approximate # of data points estimated from the plot 
% (note that because of sqrt(log(nsamp)), the result is not very sensitive to the number ; assuming ~30 points for each time point
gcm=0.5*(gcu+gcl) ; % average deviation
sc=1/(0.9*sqrt(2*log(30))) ; % scale down by approximately 0.43, i.e. roughly half of the min/max error bars
% fit
[s,delta,err]=fnrfit(t(1:end), gc(1:end), time(1:end-1), scale*dpldt',1./(gcu+gcl),1,1,1,-8.39,0.5,1)

errorbar(t+delta,gc/s,sc*gcm/s,sc*gcm/s,'k*')

ylabel('\it PC production (a.u.)', 'fontsize', 14) ;
xlabel('\it t(days)', 'fontsize', 14) ;
%xlim([0 max_t]);
xlim([0 35]);
