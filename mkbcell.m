% declare seeding Bcells :
if (~exist('qdrive') || ( exist('qdrive') && qdrive==0 ) )
 nbcr=2; % default number of clone types, unless set elsewhere (e.g. in an external "driver" script)
end
%====== conventionally, model the first AB as a stem AB (e.g. with reduced valency)
i=1;
bcr(i).affinity_min=kp93_affinity_min;
bcr(i).affinity_max=kp93_affinity_max;
bcr(i).nclass=nclass0;
if (exist(['k',num2str(i),'2'])) % check if we attempted to set 2nd affinity const. elsewhere
 eval(['bcr(i).k2=k',num2str(i),'2;']); % script to resolve at runtime
else
 bcr(i).k2=0 ; % to set here if running w/o driver
end
bcr(i).mu0=mu0;
bcr(i).lambda_scale = 1;
% set initial condition ; different ways to do this :
% simplest : populate a bcr class, e.g. xc(istep,10)=2*minimum_xi ;
% KP93 approach was to populate a single affinity class (the one with nondimensional affinity=7.5)
if (exist(['ibc',num2str(i)])) % check for external ic
 eval(['bcr(i).ixc=ibc',num2str(i),';']);
else % default
 bcr(i).ixc = zeros(1,bcr(i).nclass) ;
 bcr(i).ixc(8) = 128 * minimum_xi ; % KP93 notation
end
%=========================================
% add bcrs below as needed until ibcr=nbcr
%=========================================
while (i<nbcr ) % additional B-cells
 i=i+1;
 bcr(i).affinity_min=kp93_affinity_min; % set affinity limits
 bcr(i).affinity_max=kp93_affinity_max;
 bcr(i).nclass=nclass0;  % can set affinity discretization separately for each bcr
 if (exist(['k',num2str(i),'2'])) % for external spec
  eval(['bcr(i).k2=k',num2str(i),'2;']);
 else % default set BELOW:
  bcr(i).k2=1 ;% to model binding of two Fabs, second binding constant proportional to the first (affinity); 0 means no second fab
 end
 bcr(i).mu0 = mu0 ;% constant prefactor of mutation rate
 bcr(i).lambda_scale = 1 ;% a scaling constant for lambda that changes the ease of making advantageous mutations relative to KP93
% initial condition :
% check for external init :
 if (exist(['ibc',num2str(i)])) % check for external ic
  eval(['bcr(i).ixc=ibc',num2str(i),';']);
 else
  bcr(i).ixc = zeros(1,bcr(i).nclass) ;
  bcr(i).ixc=1*circshift(bcr(1).ixc,0,2) ; % shift to the right by arg2 places along dim arg3; concentration relative to #1
 end
end
%=======================================================================
% make sure we do not have runaway bcrs, e.g. from previous run
%
bcr=bcr(1:nbcr);
%
%========================================================================
%==== remaining initialization that should not need manual modification :
%========================================================================
for ibcr=1:nbcr
 bcr(ibcr).k2=bcr(ibcr).k2*qavid ; % nonzero only when avidity code on

% create affinity grid in log space :
 bcr(ibcr).affinity = exp ( linspace ( log(bcr(ibcr).affinity_min), log(bcr(ibcr).affinity_max), bcr(ibcr).nclass ) ) ;
% initialize arrays
 bcr(ibcr).dxc=zeros(estep,bcr(ibcr).nclass); % bcr concentration rate of change (for postpro)
 bcr(ibcr).xc=zeros(estep,bcr(ibcr).nclass); % bcr concentration arrays
 bcr(ibcr).xc(1,:) = bcr(ibcr).ixc; % apply initial condition ; NOTE that setting at index 1, which may not be general enough
 bcr(ibcr).hc=zeros(estep,bcr(ibcr).nclass); % activation (bound fraction) levels (h in KP93)
 if (qtcell)
  bcr(ibcr).hct=zeros(estep,bcr(ibcr).nclass) ;% T-cell activation model
 end
 if (qplasma)
  bcr(ibcr).xpl=zeros(estep, bcr(ibcr).nclass);
  bcr(ibcr).xmpl=zeros(estep, bcr(ibcr).nclass); % plasma cells derived from memory
  if (qab)
   bcr(ibcr).xab=zeros(estep, bcr(ibcr).nclass);
   bcr(ibcr).xmab=zeros(estep, bcr(ibcr).nclass); % ABs from memory-derived plasma cells
  end
 end
 if (qmkmem)
  bcr(ibcr).xmem=zeros(estep, bcr(ibcr).nclass);
 end
% mutation probability factors mu, lambda :
 bcr(ibcr).mu = bcr(ibcr).mu0 * ones(1,estep) ;
 bcr(ibcr).lambda = bcr(ibcr).lambda_scale* ... 
  (bcr(ibcr).affinity_max/bcr(ibcr).affinity_min)^(log(kp93_lambda)/log(7.5)/(bcr(ibcr).nclass-1));% log-space interpolation
% 1/21 : essentially lambda characterizes the distribution of AB affinities; more specifically, in KP93,
% as AB affinity goes up one class, corresponding to the increase of x7.5, the probability of finding such an AB
% goes down by 1/lambda0 (1/30); so 30 per 7.5 ; we take logs of both, interpolate linearly to "nclass" classes
% i.e. log(lambda) = log(lambda0) / log(7.5) * log(affinity(2)/affinity(1)
% use getm as a subroutine to compute transition/mutation rate matrix for this bcr ; need to set some key vars, first
 nclass=bcr(ibcr).nclass ;
 mu=bcr(ibcr).mu ;
 lambda=bcr(ibcr).lambda ;
 affinity=bcr(ibcr).affinity ;
 getm; % compute transition matrix & rate arrays for this bcr
 bcr(ibcr).mtrans=mtrans ; % store matrix
 bcr(ibcr).o=o; % also store rate arrays :
 bcr(ibcr).w=w;
 bcr(ibcr).e=e;
end
