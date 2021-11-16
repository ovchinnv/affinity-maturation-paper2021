qavid = 1; % need this to be either 1 or 0 (avidity flag) ;
qhprof = 1 ;% 1 = activated proliferation (KP93) ; 
% param below are reasonable starting poits for non-activated model: 
q=1.5;
k_prof = q * 0.85 * (1+qavid); % /day ; rate of B cell proliferation at maximum stimulation ;
k_death= q * 3 / (1+qavid); % /day ; rate of B cell death ;
%
prob_lethal=0.5 ; % mutation probability that a mutation we make is lethal ; used in getm.m
minimum_xi=1e-4 ;% concentration that corresponds to one cell (as in KP93)
bcr_per_cell=1e5 ;% number of B cell receptors per B cell
qminxi=0 ;% whether to zero out very small concentrations as in KP93 ; affects other params, not used in our paper
% affinity classes (KP93) ; we do not model explicit ABs, only how many we have with a particular affinity
% there are 8 classes in the paper ; they were indexed -3 : 4
% instead, I set affinity bounds/range, then discretize the affinity range into a prescribed # classes
kp93_affinity_min = 7.5 * 7.5 ^ (-3) ; % classes from KP93
kp93_affinity_max = 7.5 * 7.5 ^ (4 ) ;
nclass0 = 20 ; % Number of affinity classes here
adjacent_transitions=1 ; % if true, M is zero for transitions between non-adjacent states (no sig. diff. for large lambda, below)
rate_arrays=0 ; % if true, use rate arrays when adjacent_transitions = 1, rather matrix computations
% the reason including rate_arrays is that it corresponds to a sequential, non vector, fortran-like code that is a bit clearer to read
% should not materially change results here
% ==========================================
% ========= time integration parameters ==== :
it=0; % initial time, days
max_t=35; % final time, days
dt=0.01; % time step, days
t_out=1 ; % time interval between outputs
%
nstep = ceil ( ( max_t - it ) / dt ) ; % number of integration steps required to reach specified time
istep = 1 ;% current time step
estep = nstep + istep ;

% distribution parameters for computing transition matrix :
% factor by which number of variable-region sequences shrinks as affinity class increases
kp93_lambda=30 ;  % =30 to match perelson98
lambda0 = (kp93_affinity_max/kp93_affinity_min)^(log(kp93_lambda)/log(7.5)/(nclass0-1)) % log-space interpolation
% 1/21 : lambda basically characterizes the distribution of AB affinities; more specifically, in Perelson93,
% as AB affinity goes up one class, corresponding to the increase of x7.5, the probability of finding such an AB
% goes down by 1/lambda0 (1/30); so 30 per 7.5 ; we take logs of both, interpolate linearly to "nclass" classes
% i.e. log(lambda) = log(lambda0) / log(7.5) * log(affinity(2)/affinity(1)
% note that in KP93 was optimized as an oscillator but with 8 classes was about 0.1
mu0 = 0.1 ;
%============== the basic/original KP93 is above =============
% additions below
%==============
% maximum GC size
gcmax=inf ;% no restriction on total size
clonemax=5e3 * minimum_xi ; % maximum clone size ; equivalent to gcmax for single-AG GCs
% decay rate of IC bound antigen :
% in this simplest model, there is only one decay parameter, no growth parameter
kalpha_death0 = 0 ;% to turn off AG decay
kalpha_death0 = 4.63e-4 * 24 ; % Rundell 1998
%
qplasma=1 ; % whether long lived plasma cells (pl) are modeled
qab=1;% whether antibodies are modeled
qmkmem=1 ;% whether to create memory cells
qtcell=0 ; % T cell help model ; 0 means that T help is simply the fraction bound to antigen ; 1 means model T-cell concentration
if (~exist('hcexp'))
 hcexp=0.7; %exponent of activation function, which is hc or hct, depending on qtcell (1 is the default that corresponds to KP93)
end
%
if (~exist('hcscale'))
 hcscale=0.94 ; % scaling for the activation function h
end
% note that we assume that the range of h above should be 0-1, otherwise will get unexpected results ( high B-cell growth )
% GC exit activation function : C x [ h^p (1-h)^q ] ; note that an extra factor of (1-h) is added in the integrator integ2.m
% note that it assumes that the range of 1 is [0-1] ; this basically a splitting function, which thus needs to be normalized
if (qplasma || qmkmem)
 pexit=1.; qexit=0 ; pqnorm=(pexit+qexit)^(pexit+qexit)/(pexit^pexit*qexit^qexit) ; % computationally "bad" for large p,q, BUT stable near 0
 if (qplasma)
  Cpl=0.7 ;% addtl constant in differentiation rate to plasma cells
  kpl_death=0.0336 ; % per day ; from R98
% no current memory source
% antibody production by plasma cells
  if (qab)
   kab_death = log(2) / 10  / bcr_per_cell ; % per day (R98) ; Zhang13 have log(2)/10 = 0.069 / day
   kab_prod = 35e3  / bcr_per_cell; % estimated in Methods of the paper
  end
 else % qplasma
  Cpl=0; % make sure this is set
  qab=0; % need plasma cells to make Abs
 end
% memory
 if (qmkmem)
  Cmem=0.3 ;% addtl constant in differentiation rate to memory cells
  kmem_death=0.02 ; % per day ; from R98 : ~ 8.5e-4/h x 24 h/d  ; note how much smaller this is than Bcell diff rate
  if(qplasma)
   kmpl_prod=0.18  % 0.17 - 0.18 gives a reasonable (ad hoc) fit to Weisel et al '16 data
   Cmemab=1; % % ABs derived from memory cells compete with this scaling (0 turns off)
  end
 else
  Cmem=0;
 end
%
else
 Cpl=0;Cmem=0;
 qab=0; % need plasma cells to make abs
end
%
%===================================================
% == call Bcell spec :
mkbcell ;
% == call Antigen spec :
mkag ;
% == T cell model spec :
mktcell ;
%
% make sure we remove unneeded vars : mtrans, o, w, e to avoid using them mistakenly (params are now stored in bcr struc)
%
clear mtrans
clear o
clear w
clear e
%
