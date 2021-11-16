% whether to initialize (otherwise, continue from prev. state) :
if (~exist('qinit'))
 qinit=1;
end
%
if (qinit)
% errcalc='driver;ifig=1;gcexp;rerr1=err/sum(gc.^2);gcmbc;rerr2=err/sum(gc.^2);err=rerr1+3*rerr2 ;'% script to run to compute error
 errcalc='driver;ifig=1;gcexp;rerr1=err/sum(gc.^2);gcmbc;rerr2=err/sum(gc.^2);gcplc;rerr3=err/sum(gc.^2);err=rerr1+3*rerr2+2*rerr3;'% script to run to compute error
 clear rave;
 clear opt;
 maxoptiter=100;
 i=0; % opt var index
% create optimization structure
% specify structure
 i=i+1;
 opt(i).name='hcexp'; % must match parameter name in error evaluator
 opt(i).ival=1 ; % initial value
 opt(i).dh=0.01 ;% FD step in computing derivative
 opt(i).step=0.02 ;% evolution step in gradient descent
 opt(i).maxchange=0.1 ;% maximum change between iterations
 opt(i).minval=0.5 ; % minimum allowed value (-inf to disable)
 opt(i).maxval=1. ; % maximum allowed value (inf to disable)
%
 i=i+1;
 opt(i).name='hcscale';
 opt(i).ival=1 ;
 opt(i).dh=0.01 ;
 opt(i).step=0.02 ;
 opt(i).maxchange=0.1 ;
 opt(i).minval=0.5 ;
 opt(i).maxval=1. ;
%% generic inits below
 for j=1:i
  opt(j).val=opt(j).ival ;
 end
 npar=length(opt); % number of parameters to optimize
 xpar=zeros(0,npar);
 rr=zeros(1,0);
%
 qinit=0;
end
%
fmt='%.16f'; % format for num2str
nrepe = 1 ; % number of repetitions to compute error (i.e., if stochastic)

for optiter=1:maxoptiter
 figure(1) ; clf ; figure(2) ; clf;
 fprintf('================== Iteration %d\n', optiter);
% set initial values to make sure that all perems are defined:
 for ipar=1:npar
  eval( [opt(ipar).name,'=',sprintf(fmt,opt(ipar).val),';']);
 end
% compute error derivative by FD
 for ipar=1:npar
  dd=[-opt(ipar).dh +opt(ipar).dh];
%
  for k=1:2 % compute diff
   eval( [opt(ipar).name,'=',sprintf(fmt, opt(ipar).val + dd(k) ),';']);
% custom:
   for irepe=1:nrepe
    eval(errcalc); % run err calc
%    return
    rerrs(k, irepe)=err
   end
   rerr=mean(rerrs,2);
% custom
  end
  drerr(ipar)= (rerr(2)-rerr(1))/(dd(2)-dd(1))
  rave(ipar) = 0.5*(rerr(1)+rerr(2))
 end % over all params
%
% compute new parameter values
%
 for ipar=1:npar
  dh = min( opt(ipar).maxchange, opt(ipar).step * abs(drerr(ipar))); % evolution size
  opt(ipar).val = opt(ipar).val - dh*sign(drerr(ipar));
  opt(ipar).val = min(opt(ipar).maxval, max( opt(ipar).val, opt(ipar).minval )); % make sure to keep inside a prescribed interval
 end

 xpar=[xpar ; [opt(:).val] ];

 rr=[rr mean(rave)]; % keep track of error

 fprintf('%s %d\n', ' === iteration',optiter);
 fprintf( ['%s ', repmat('%17.12f', [1,npar]), '\n'], 'New parameters: ',[opt(:).val]);
%
end % iterations
%
savename=['fit-',date,'.mat'];
save(savename, '-mat')
