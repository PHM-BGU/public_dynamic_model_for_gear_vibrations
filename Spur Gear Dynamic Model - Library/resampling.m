function xNew = resampling(xOld, NPtsNew, interpMethod)
%{
% Description:
% This function operates on a vector/matrix and performs
% resampling based on the specified number of points.
% The default interpolation method is set to 'spline'. However, if
% 'xOld' contains sharp discontinuities, such as variations in tooth
% width in the case of a tooth face fault, 'PCHIP' would be a preferable 
% choice to avoid generating artifacts around these sharp edges.
% =====
% Inputs:
% * xOld - the original vector.
% * NPtsNew - the number of points in the new resolution.
% * interpMethod - optional, interpolation method (spline/ PCHIP/ linear).
% =====
% Outputs:
% * xNew - the resampled vector in new resolution.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Check for interpolation method %%
if nargin == 2
    interpMethod = 'spline' ;
end

%% Resampling by interpolation %%
Npt_old = length(xOld) ;
coarseVctr = linspace(0, 1, Npt_old)' ; % neutral vctr for interpolation (coarse)
fineVctr = linspace(0, 1, NPtsNew)' ; % neutral vctr for interpolation (fine)
xNew = interp1(coarseVctr, xOld, fineVctr, interpMethod) ;

end % of function 'resampling'