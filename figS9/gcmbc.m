% compute mbc production rate from simulation data
% from integ2 :
% hex=pqnorm * (hactivation.^pexit) .* ((1 - hactivation).^qexit); % cell exiting proliferation/mutation cycle
% dexdt = (1-hactivation) .* hex .* xcorr ;
% bcr(ibcr).xmem(istep+1,:) = (1 - dt*kmem_death) * bcr(ibcr).xmem(istep,:) + dt * Cmem * dexdt ;
colors={'b', 'g', 'r', 'm', 'c', 'k'};
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
scale=1/(75) ;
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
gcu=d(2:3:end,2)-gc; % max value
gcl=gc-d(3:3:end,2); % min value
toff = 1;
sc=0.5; % crude scale of minimum & maximum values to get estimate of standard deviation

%[s,delta,err]=fnrfit(t(1:end), gc(1:end), time(1:end-1), scale*dmemdt',1,1,0,1,0,0.5,1)
%[s,delta,err]=fnrfit(t(1:end), gc(1:end), time(1:end-1), scale*dmemdt',1./(gcu+gcl),1,0,1,1.5,0.5,1)
[s,delta,err]=fnrfit(t(1:end), gc(1:end), time(1:end-1), scale*dmemdt',1./(gcu+gcl),1,0,1,1,0.5,1)
errorbar(t+delta,gc/s,sc*gcl/s,sc*gcu/s,'k*')

ylabel('\it MBC production (a.u.)', 'fontsize', 14) ;
xlabel('\it t(days)', 'fontsize', 14) ;
xlim([0 35]);
