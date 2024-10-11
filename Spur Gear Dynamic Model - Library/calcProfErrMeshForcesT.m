function profErrMeshForcesT = calcProfErrMeshForcesT(errsMtx, KeqStruct, interpMethod, dt, rpsSig)
%{
% Description:
% This function computes mesh forces induced by profile errors, such as
% surface quality imperfections or fault displacements originating from
% involute distortions. The input comprises a vecotr or a matrix of profile
% errors errsMtx, where each column corresponds to a distinct meshing pair, and a
% structure KeqStruct representing the equivalent stiffness of a healthy tooth pair.
% Optionally, KeqStruct can also provide the stiffness of faulty tooth pairs.
% First, the profile errors are interpolated to match the finer resolution
% of the equivalent stiffness. Subsequently, the error for each tooth pair
% is scaled by its corresponding stiffness. The resulting meshing forces
% for all tooth pairs are concatenated iteratively, ensuring proper
% overlap based on the contact ratio (ε). In the final step, the force
% signal in the cycle domain is converted to the time domain.
% =====
% Inputs:
% * errsMtx - [mm], profile errors, can be either a vector or a matrix.
%   If a vector is provided, it is assumed that only the first tooth pair
%   is associated with this vector, and the mesh forces induced by the
%   remaining tooth pairs are considered zero. In the case of a matrix,
%   each column is associated with a distinct tooth pair.
% * KeqStruct - a strucT containing equivalent stiffness & essential information:
%   # KeqStruct.KeqHealthy - [N/mm], the equivalent mesh stiffness of a
%   healthy tooth pair (lenK x 1).
%   # KeqStruct.KeqFaulty - [N/mm], the equivalent mesh stiffness of
%   the first M faulty tooth pairs (lenK x M) (optional).
%   # KeqStruct.gmCycDur - [cycle], the cycle duration of a complete
%   gearmesh period (1/z, usually zIn).
%   # KeqStruct.epsilon - Contact ratio (ε),indicating the number of
%   gearmesh cycles elapsed from engagement to separation. Keq elapses ε
%   gearmesh cycles.
%   # KeqStruct.dCyc - [cyc], resolution of Keq in cycle, which is also the
%   desired resolution of errsMtx in cycle after interpolation.
%   # KeqStruct.z2Repeat - number of repetitions of Keq required to create
%   a matrix of stiffnesses corresponding for each pair.
% * interpMethod - the interpolation method for errsMtx (spline/PCHIP, etc).
%   Sharp signals require 'PCHIP', while smooter signals 'spline'.
% * [rpsSig, dt] - input speed signal [Hz] and the time resolution [s].
% =====
% Output:
% * profErrMeshForcesT - [N], The mesh forces signal in time domain.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Extract all the relevant information from KeqStruct %%
[KeqHealthy, gmCycDur, epsilon, dCyc, z2Repeat] = deal( ...
    KeqStruct.KeqHealthy, KeqStruct.gmCycDur, KeqStruct.epsilon, ...
    KeqStruct.dCyc, KeqStruct.z2Repeat) ;

%% Construct a matrix of equivalent stiffness for each tooth pair %%
KeqMtx = repmat(KeqHealthy, 1, z2Repeat) ;
try KeqMtx(:, 1:size(KeqStruct.KeqFaulty, 2)) = KeqStruct.KeqFaulty ;
end % of try

%% Interpolate the errors matrix to match the KeqMtx resolution %%
[lenErrsSig, NumToothPairs] = size(errsMtx) ;
cycPtsCoarse = linspace(0, gmCycDur*epsilon, lenErrsSig)' ; % for ε cycles
cycPtsFine = (0 : dCyc : cycPtsCoarse(end))' ;
errsMtxFine = interp1(cycPtsCoarse, errsMtx, cycPtsFine, interpMethod) ;

%% If errsMtx is a vector, only the first tooth pair contributes %%
if NumToothPairs == 1
    gmCycNum = z2Repeat ;
else
    gmCycNum = NumToothPairs ;
end % of if

%% Calculate profile errors-induced mesh forces iteratively in cycle domain %%
numSampPerCycle = ceil(gmCycDur/dCyc) ; % number of points in one gearmesh cycle
profErrMeshForcesCyc = zeros(numSampPerCycle*gmCycNum, 1) ; % pre-allocation
for ii = 1:NumToothPairs %% utilizing modulo to capture cyclic nature
    profErrMeshForceCurrPair = KeqMtx(:, mod(ii-1, z2Repeat)+1) .* errsMtxFine(:, ii) ;
    startInd = 1 + (ii-1)*numSampPerCycle ;
    endInd = startInd + length(profErrMeshForceCurrPair) - 1 ;
    currPairInds = mod((startInd:endInd)-1, length(profErrMeshForcesCyc)) + 1 ;
    profErrMeshForcesCyc(currPairInds) = profErrMeshForcesCyc(currPairInds) +  profErrMeshForceCurrPair ;
end % of for

%% Convert the profile error mesh forces signal from cycle domain to time domain %%
profErrMeshForcesT = convertCyc2Time(dCyc, profErrMeshForcesCyc, dt, rpsSig) ;

end % of function 'calcProfErrMeshForcesT'