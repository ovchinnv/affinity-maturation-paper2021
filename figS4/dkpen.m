% make a plot that shows dk penalty of ste vs others ABs
%figure(1, 'position', [100 100 350 500])
%clf

nrun=10000;

% distribution parameters :
s=1;
mu=2;

s=sqrt(0.5); mu=1. ;

dlk=(s*randn(1,nrun)+mu)  % sample from gaussian in log space centered on 10^mu

[n,x]=hist(dlk,25) ;
bar(x,n/sum(n)/(x(2)-x(1)),'facecolor','white')


xlim([mu-2 mu+2])

%set(gca, 'xscale', 'log')

% plot gaussian pdf : 
xx=linspace(-5,5,100);
pdf=exp(-0.5*((xx-mu)/s).^2)/sqrt(2*pi*s*s);
hold on ; 
plot(xx,pdf,'k', 'linewidth',1)
set(gca, 'fontsize', 11)
l=xlabel( 'K_{eq}^{11} / K_{eq}^{i >1,1}', 'fontsize', 20 )
set(l, 'fontsize',14)

xl=get(gca,'xticklabel')
clear xlnew;
for i=1:length(xl)
 xlnew(i)=str2num(char(xl(i)));
end
xlnew=10.^xlnew / 1000;
xlnew=1./xlnew(end:-1:1) ; % invert
set(gca, 'xticklabel', xlnew )

