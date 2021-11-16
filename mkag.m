%============================================
% declare antigens
nag=nbcr; % could be arbitrary, but needs to be set in coordination with #bcrs (i.e. the binding partners)
% initialize AG occlusion matrix
occlusion=eye(nag); % initialize to independent antigens
%
% could specify manually, or in a loop, as below ( identical ags )
for iag=1:nag
 ags(iag).kalpha_death=kalpha_death0; % set to a global constant
 ags(iag).bcr=iag; % index into bcr array, for bcr(s) that binds to this ag ; note that in the our model, a bcr only binds to 1 ag
%
 if (exist(['agc',num2str(iag)])) % for external spec
  eval(['ags(iag).iagtot=agc',num2str(iag),';']);
 else % default set below:
  ags(iag).iagtot=1 ; % initial total antigen concentration
 end
end
%
% specify simple occlusion matrix
if (exist('occl'))
 oc=occl;
else
 oc=0.9 ; % constant
end
for iag=1:nag
 for jag=iag+1:nag
  occlusion(iag,jag)=oc ; % 1 -- maximum occlusion
  occlusion(jag,iag)=oc ;
 end
end
%
% make sure we do not have runaway ags, eg, from previous calc.
%
ags=ags(1:nag);
%========================================================================
%==== remaining initialization that should not need manual modification:
%========================================================================
for iag=1:nag
 ags(iag).alphatot = zeros(1,estep) ; % total antigen
 ags(iag).alphatot(istep) = ags(iag).iagtot ;
 ags(iag).alpha = ags(iag).alphatot ;   % free antigen (intial guesses)
end

%== normalize occlusion matrix :
for iag=1:nag
 for jag=1:nag
  occlusion(iag,jag)=occlusion(iag,jag) * min(1.0, ags(iag).iagtot / ags(jag).iagtot ); % normalization to avoid negative [a]
 end
end
%
