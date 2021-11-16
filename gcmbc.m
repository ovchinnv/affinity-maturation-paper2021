% compute mbc production rate from simulation data
% from integ2 :
% hex=pqnorm * (hactivation.^pexit) .* ((1 - hactivation).^qexit); % cell exiting proliferation/mutation cycle
% dexdt = (1-hactivation) .* hex .* xcorr ;
% bcr(ibcr).xmem(istep+1,:) = (1 - dt*kmem_death) * bcr(ibcr).xmem(istep,:) + dt * Cmem * dexdt ;
colors={'k', 'g', 'r', 'm', 'c', 'b'};
if (~exist('qmkmem'))
 qmkmem=0;
end
if (qmkmem==0)
 return
end
%
if (~exist('ifig')) ; ifig=1 ; end
%
ifig=ifig+1;figure(ifig)
leg=[];
hold on ; box on ;
scale=1/(75) ; % arbitrary scaling (exp data in arbitrary units)
% to use finite difference with existing MBC counts:
time = dt * [0:nstep];
for iag=1:nag
 for ibcr=ags(iag).bcr
   tmem=sum(bcr(ibcr).xmem,2);
   difftmem=sum(diff(bcr(ibcr).xmem),2);
   dmemdt=(difftmem/dt + kmem_death * tmem(1:end-1))/minimum_xi ;
   plot(time(1:end-1),scale*dmemdt,[char(colors(ibcr)),'-'] ) ;
   leg=[leg { ['AG #',num2str(iag),'; MBC #',num2str(ibcr)] } ];
 end
end
% plot exp weisel 16 data from pelissier 20
d=load('../data/pelissier4e.dat');

t=d(1:3:end,1);
gc=d(1:3:end,2); % average
gcu=d(2:3:end,2)-gc; % maximum upper deviation
gcl=gc-d(3:3:end,2); % maximum lower deviation
% crude scale of minimum & maximum values to get estimate of standard deviation :
% to get the SD estimate (sig), assume that the data are Gaussian distributed and use the estimate :
% x_max ~ 0.9 * sig * sqrt( 2 * log (nsamp) ) ; where nsamp is the approximate # of data points estimated from the plot 
% (note that because of sqrt(log(nsamp)), the result is not very sensitive to the number ; assuming ~30 points for each time point
% this is a standard relation obtained from Jensen's inequality and the definition of a moment generating function ;
% the x0.9 is an empirical correction for smaller samples
gcm=0.5*(gcu+gcl) ; % average deviation
sc=1/(0.9*sqrt(2*log(30))) ; % approximately 0.43, i.e. roughly half of the min/max error bars
%
%[s,delta,err]=fnrfit(t(1:end), gc(1:end), time(1:end-1), scale*dmemdt',1,1,0,1,0,0.5,1)
%[s,delta,err]=fnrfit(t(1:end), gc(1:end), time(1:end-1), scale*dmemdt',1./(gcu+gcl),1,0,1,1.5,0.5,1)
[s,delta,err]=fnrfit(t(1:end), gc(1:end), time(1:end-1), scale*dmemdt',1./(gcu+gcl),1,0,1,1,0.5,1)
errorbar(t+delta,gc/s,sc*gcm/s,sc*gcm/s,'k*')

ylabel('\it MBC production (a.u.)', 'fontsize', 14) ;
xlabel('\it t(days)', 'fontsize', 14) ;
xlim([0 35]);
