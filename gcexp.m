% plot wittenbrink 11 data from pelissier 20
colors={'k', 'g', 'r', 'm', 'c', 'b'};
d=load('../data/wittenbrinks1b1.dat');

t=d(1:3:end,1);
gc=d(1:3:end,2);
gcu=d(2:3:end,2)-gc;
gcl=gc-d(3:3:end,2);
toff = 0;
%
if (~exist('ifig')) ; ifig=1 ; end
%
figure(ifig) ; hold on ;
time = dt * [0:nstep];
for iag=1:nag
 for ibcr=ags(iag).bcr
  tx=sum(bcr(ibcr).xc,2)/minimum_xi ;
  plot(time,tx,[char(colors(ibcr))] ) ;
#  leg=[leg { ['AG #',num2str(iag),'; BCR #',num2str(ibcr)] } ];
 end
end

%[s,delta,err]=fnrfit(t, gc, time(:), tx(:), 0, 1, 1, -3, 0.5)
%[s,delta,err]=fnrfit(t, gc, time(:), tx(:), 1, 0, 1, 1, -0.5, 0.5)
[s,delta,err]=fnrfit(t, gc, time(:), tx(:), 1./(gcu+gcl), 0, 1, 1, -0.5, 0.5)

errorbar(t+delta,gc/s,gcl/s,gcu/s,'k*')

