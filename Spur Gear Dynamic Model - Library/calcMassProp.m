function wheelInfo = calcMassProp(wheelInfo)
%{
% Description:
% This function calculates the mass and tensor of inertia of a gear wheel
% as if it was a rigid cylindrical body with a radius of the pitch circle
% and the height of the tooth width. The function takes into account the 
% option of pre-defined mass properties ('massProps') provided by the user.
% =====
% Inputs:
% * wheelInfo - relevant field in this structure:
%   # wheelInfo.module - [mm].
%   # wheelInfo.z - number of teeth.
%   # wheelInfo.W - [mm], tooth width.
%   # wheelInfo.rho - [kg/m^3], density.
% =====
% Outputs:
% * wheelInfo - relevant field in this structure:
%   # wheelInfo.m - [kg], mass.
%   # wheelInfo.rotI - [kg m^2] - rotational moment of inertia.
%   # wheelInfo.polarJ - [kg m^2] - polar moment of inertia.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% step 1: define the pitch radius and tooth width %%
R = 0.5 * wheelInfo.module * wheelInfo.z * 1e-3 ; % from [mm] to [m], pitch radius.
W = wheelInfo.W * 1e-3 ; % from [mm] to [m], tooth width

%% step 2: calculate the mass or allow the user to set it %%
name = wheelInfo.wheel ;
try
    wheelInfo.m = wheelInfo.massProps.(name).m ;
catch
    wheelInfo.m = wheelInfo.rho * pi * R^2 * W ;
end % of try-catch

%% step 3: calculate the tensor of inertia or allow the user to set it %%
try
    wheelInfo.rotI = wheelInfo.massProps.(name).rotI ;
    wheelInfo.polarJ = wheelInfo.massProps.(name).polarJ ;
catch
    wheelInfo.rotI = wheelInfo.m * (3*R^2 + W^2) ;
    wheelInfo.polarJ = 0.5 * wheelInfo.m * R^2 ;
end % of try-catch

%% step 4: remove the mass properties field (optional) %%
try wheelInfo = rmfield(wheelInfo, 'massProps') ;
end % of try

end % of function 'calcMassProp'