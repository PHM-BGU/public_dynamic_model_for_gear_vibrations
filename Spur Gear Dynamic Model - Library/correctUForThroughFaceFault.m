function U = correctUForThroughFaceFault(U, faultInds, F1, F2)
%{
% Description:
% This function corrects the strain energy U for the case of ThroughFaceFault.
% The indicies in U related to the fault (where there is no actual contact)
% are overwritten with the value of the proper fault edge, and then retaining
% the original α-dependent term F1F2(α) by division and multiplication.
% According to ThroughFaceFault, fault indices are divided into two ranges:
% 1. initiation-transition: the undamaged wheel rotates along the starting
%    edge of the damaged wheel until the middle point.
% 2. transition-exit: the undamaged wheel rotates along the ending edge of
%    the damaged wheel until it exits the defective region.
% =====
% Inputs:
% * U - Potential energy vector.
% * faultInds - structure with indices corresponding to fault's start, middle, and end.
% * [F1, F2] - Meshing force componenets (either axial and/or shear).
% =====
% Outputs:
% * U - The modified U vector according to the ThroughFaceFault modeling.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
%% step 1: input check %%
if nargin == 3
    F2 = F1 ;
end % of if

%% step 2: load the fault indices and F1F2(α) %%
F1F2 = F1.*F2 ;
if faultInds.end > faultInds.start
    initTransInds = faultInds.start:faultInds.mid ;
    transExitInds = (faultInds.mid+1):faultInds.end ;
else
    initTransInds = faultInds.mid:faultInds.start ;
    transExitInds = faultInds.end:(faultInds.mid-1) ;
end % of if

edgeInds = [faultInds.start, faultInds.end] ;
indRanges = {...
    initTransInds ; ...  % initiation-transition
    transExitInds } ;    % transition-exit

%% step 3: overwrite U_fault with U_edge and retain F1F2(α_fault) %%
for ii = 1:length(edgeInds)
    U(indRanges{ii}) = ...
        U(edgeInds(ii)) .* F1F2(indRanges{ii}) ./ F1F2(edgeInds(ii)) ;
end % of for ii
end % of function 'correctUForThroughFaceFault'