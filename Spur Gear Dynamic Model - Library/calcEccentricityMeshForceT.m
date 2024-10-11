function eccMeshForcesT = calcEccentricityMeshForceT(eccentricityLevel, gmsStruct, tr, alpha, dt, rpsSig)
%{
% Description:
% This function calculates the mesh forces induced by eccentricity.
% The calculated eccentricity takes both wheels (in and out) into account.
% Meshing force is computed as the product of eccentricity and gear mesh
% stiffness, both over time.
% =====
% Inputs:
% * eccentricityLevel - [mm], eccentricity level.
% * gmsStruct.gmsCyc - [N/mm], gearmesh stiffness in cycle.
% * gmsStruct.dCycleCoarse - [cyc], resolution of the gms in cycle.
% * tr - transmission ratio (zOut/zIn).
% * alpha - [rad], nominal pressure angle.
% * [rpsSig, dt] - In speed signal [Hz] and time resolution [s].
% =====
% Outputs:
% * eccMeshForcesT - [N], The mesh forces vector in time.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
%% Input Check %%
if eccentricityLevel == 0
    eccMeshForcesT = 0 ;
    return
end % of if

%% Calculate the eccentricity signal in time %%
cycMotor = 2*pi*cumsum(rpsSig*dt) ;
modulationSigIn = sin(alpha + cycMotor) ;
modulationSigOut = sin(alpha + cycMotor/tr) ;
eccentricityT = eccentricityLevel .* (modulationSigOut - modulationSigIn) ;
%% Calculate the gearmesh stiffness in time %%
gmsT = convertCyc2Time(gmsStruct.dCycle, gmsStruct.gmsCyc, dt, rpsSig) ;
%% Calculate the mesh forces %%
eccMeshForcesT = gmsT .* eccentricityT ;
end % of function 'calcEccentricityMeshForceT'