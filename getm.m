% compute mutation matrix (eqs. 17 -- 19 in KP93)
% note that we can compute all matrix entries here ; i.e. even though per KP93, only adjacents are large enough to be significant

mtrans=zeros(nclass); % nclass x nclass matrix of affinity transitions (like a Markov transition matrix)
% use exponential (PK93)
for i=1:nclass
  for j=i+1:nclass - max(0,(nclass-(i+1))) * (adjacent_transitions)
    dclass=j-i;
% pre-compute using logs
    mlog = dclass * log( mu(istep) * (1-prob_lethal) ) - sum ( log (1:dclass ) ) - mu(istep);
%
    mtrans(i,j) = exp(mlog) ;
    mtrans(j,i) = mtrans(i,j) / (1+lambda^(-2*dclass)); % must come before next line, otherwise overwriting m_ij !
    mtrans(i,j) = mtrans(i,j) / (1+lambda^(2*dclass));
  end
%  mtrans(i,i) = exp(-mu(istep)); % from PK93 paper, but does not conserve probability discretely (implying a lower division rate)
  mtrans(i,i)=1-sum(mtrans(i,:)) ; % do I need the 1 ? Yes, because sum(dx_i/dt) = kp * sum(x_i/dt), which expresses the fact
% that all of the cells are growing at the specified rate
end
% populate rate arrays ; only needed for adjacent transitions (see param.m & integ2.m)
o(1:nclass)=diag(mtrans)-1; % see integrator for why the difference
w(2:nclass)=diag(mtrans,1);
e(1:nclass-1)=diag(mtrans,-1);
