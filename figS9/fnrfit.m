function [s,delta,err]=fnrfit(texp, fexp, t, f, w=1, qscale=1, qshift=1, s0=1, delta0=0, dtgapscale=0.25, qscalegt0=0)
% NR iterative best fit of fexp(text) to f(t); fexp can be optionaly scaled an texp can be shifted
% parameters :
% w : weights for squared error
% qscale : whether to allow scaling of fext to improve fit
% qshift : whether to allow shifting of texp to improve fit
% s0 : initial scaling guess
% delta0 : initial shift guess
% dtgapscale : shifting can reduce overlap between t and texp ; this sets maximum gap between t and texp boundaries
% qscalegt0 : whether to disallow nonpositive scale
%
if (numel(w)==1)
 wgts=ones(size(texp));
else
 wgts=w/sum(w)*numel(texp);
end
%
if (~qscale)
 if (~s0)
  s0=1;
 end
else
 qscale=true;
end
if (~qshift)
 if (~delta0)
  delta0=0;
 end
else
 qshift=true;
end
flags=[qscale ; qshift];

% set support bounds
dtgap = dtgapscale * (texp(end)-texp(1)); % experimental support size
% iterate :
maxiter=100;
nrstep=0.25 ; % scaling factor for NR
sdstep=0.01 ; % optimization step in SD
iter=1; % iteration counter
s=s0 ;
delta=delta0;
% difficult to say which interpolation is the best; beware of distant extrpolation
method='linear';
method='spline';
%method='cubic';
%method='nearest';

qnr=1 ; % newton-raphson
qsd=1 ; % steepest decent

dt=(t(2:end)-t(1:end-1));
dtc=0.5*(dt(1:end-1) + dt(2:end));
odt=1./dt;odtc=1./dtc;

do
% compute derivatives
 df  = (f(2:end) - f(1:end-1)).*odt ;% 1st derivative
 dfc = 0.5*(df(1:end-1) + df(2:end)); % 1st derivative at the same t as f
 d2f = (df(2:end) - df(1:end-1)).*odtc;% 2nd derivative at the same t as f

 ffit = interp1(t-delta,f,texp,method,'extrap') ; % extrapolation avoids NaN when support is shifted too far away 
 dffit= interp1(t(2:end-1)-delta,dfc,texp,method,'extrap') ; %note smaller support here and below
%
 ediff = (fexp-s*ffit) ;
 err=sum(wgts.*(ediff.^2)) ;
 ediff = ediff.*wgts; % scale to avoid multiplications below
% F  = [  sum(ediff.*ffit) ; s*sum(ediff.*dffit) ];
% scale s factored out and removed in second equation (matrix not symmetric):
 F  = [ sum(ediff.*ffit) ; sum(ediff.*dffit) ];
 if (qnr)
  d2fit= interp1(t(2:end-1)-delta,d2f,texp,method,'extrap') ;
%  DF = [        -sum(wgts.*ffit.^2)      prod(flags)*sum(wgts.*dffit.*(fexp-2*s*ffit)) ; ... % prod(flags)=0 if scale or shift fixed
%     prod(flags)*sum(wgts.*dffit.*(fexp-2*s*ffit)) s*sum(d2fit.*ediff - s*wgts.*dffit.^2) ];
%  DF = [        -sum(wgts.*ffit.^2)      prod(flags)*sum(wgts.*dffit.*(fexp-2*s*ffit)) ; ... % prod(flags)=0 if scale or shift fixed
%            0                                      s*sum(d2fit.*ediff - s*wgts.*dffit.^2) ];
%  DF(2,1)=DF(1,2);
% scale s factored out and removed in second equation (note that matrix not symmetric):
   DF = [ -sum(ffit.^2)           prod(flags)*sum(dffit.*(fexp-2*s*ffit)) ; ...
          -prod(flags)*sum(wgts.*ffit.*dffit) sum(d2fit.*ediff - s*wgts.*dffit.^2) ] ;
%
 dpar = - nrstep * (DF\F).*flags ; % actually less accurate in some cases
% (DF\F).*flags
% if(any(isnan(DF)))
%  warning "DF is NaN"
%  break
% end
% if(any(isnan(F)))
%  warning "F is NaN"
%  break
% end
% (inv(DF)*F).*flags
%  dpar = - nrstep * (inv(DF)*F).*flags;  % equivalent but less stable for near singular DF
 elseif (qsd)
  dpar = + sdstep/max(1,err) * F.*flags;  % note sign, because of the way error is written
 end
 dpar ;
 s = s + dpar(1);
 if (qscalegt0)
  if(s<=0)
   warning "Negative scale encountered, resetting guess"  
   s = s - 0.99 * dpar(1) ;
   iter=iter-1;
  end
 end
 if (qshift)
  delta = delta + dpar(2) ;
% check support after delta modified
  tleft=t(1)-delta; %left boundary
  tright=t(end)-delta; %right boundary
  if (tleft > texp(1))
   if( (tleft-texp(1))>dtgap )
    warning "Left support boundary is too high, resetting guess"  
    s = s - 0.9 * dpar(1) ;
%   delta = delta - 0.99 * dpar(2) ; % adhoc
    delta=delta0 ;
%    iter=iter-1;
   else
    warning "Left support boundary is high, continue"
   end
  end
%
  if (tright < texp(end))
   if ( (texp(end)-tright)>dtgap )
    warning "Right support boundary is too low, resetting guess"  
    s = s - 0.9 * dpar(1) ;
%   delta = delta - 0.99 * dpar(2) ; % adhoc
    delta = delta0;
%    iter=iter-1;
   else
    warning "Right support boundary is low, continue"
   end
  end
 end % qshift
%
% plot(texp,s*ffit,':') ;
 iter=iter+1;
 if (iter>maxiter)
  warning "Maximum number of iterations exceeded"
  break
 end
 if (~any(flags)) % nothing to do -- both flags are zero
  break
 end
% norm(F.*flags)
until (norm(F.*flags)<1e-2)
%plot(t-delta,s*f,'k-') ;
err=sum((s*ffit-fexp).^2);
%plot(texp,s*ffit,'rx', 'linewidth',1 ) ;
