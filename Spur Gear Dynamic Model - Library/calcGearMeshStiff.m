function gmsStruct = calcGearMeshStiff(inWheelInfo, outWheelInfo, ...
    simParam, inWheelDefectInfo, outWheelDefectInfo)
%{
% Description:
% This function calculates the gear mesh stiffness (gms) as a function of
% the cycle by a proper fusion of the equivalent stiffnesses of all pairs.
% The shape of the gms is a square-wave-like curve, as illustrated below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       *******     *******     *******     %
%       *     *     *     *     *     *     %
%       *     *     *     *     *     *     %
%       *     *     *     *     *     *     %
%       *     *     *     *     *     *     %
% *******     *******     *******     ***** %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====
% Inputs:
% * (in/out)WheelInfo - a structure with all the information of the tooth.
% * (in/out)WheelDefectInfo - a structure with all the information about the fault.
% * simParam - a structure with the simulstion parameters, relevant
%   fields for this function:
%   # simParam.dCycGMSFine - the fine resolution in cycle of the gms
%   # simParam.gmsCoarseInterval - the ratio between the fine and
%   coarse resolutions of the gms.
% =====
% Outputs:
% * gmsStruct - A structure with the gms and the following information:
%   # gmsStruct.fine - [N/mm] - The gms vector with the finer resolution.
%   # gmsStruct.coarse - [N/mm] - The gms vector with the coarser resolution.
%   # gmsStruct.cycVctrCoarse - The cycle vector of the coarse gms.
%   # gmsStruct.alphaCyc - [rad], The pressure angle in cycle (coarse).
%   # gmsStruct.initContInd = The initial contact index.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
    %}
    
    %% Calculate the equivalent stiffness of a healthy tooth pair %%
    KeqHealthy = calcKeq(inWheelInfo, outWheelInfo) ;
    gmCycDur  = 1/inWheelInfo.z ; % [cyc] (~gearmesh period time)
    lenGMCyc = ceil(gmCycDur/simParam.dCycGMSFine) ; % number of points in one gearmesh cycle
    cycPtsCoarse = linspace(0, gmCycDur*inWheelInfo.epsilon, length(KeqHealthy))' ; % for ε cycles
    cycPtsFine = (0:simParam.dCycGMSFine:cycPtsCoarse(end))' ;
    KeqHealthy = interp1(cycPtsCoarse, KeqHealthy, cycPtsFine, 'spline') ;
    
    %% Calculate the gms of a healthy pair for one gearmesh cycle %%
    % Clarification: In one gearmesh cycle, the total gms is the sum of the
    % remaining stiffness from the separating tooth pair and the new
    % stiffness contributed by the engaging tooth pair.
    KeqSeparatingPair = [KeqHealthy(lenGMCyc+1:end) ; zeros(lenGMCyc,1)] ;
    KeqOneGMCyc = KeqSeparatingPair(1:lenGMCyc) + KeqHealthy(1:lenGMCyc) ;
    
    %% Repeat and concatenate Keq for a complete cycle of the shaft (z gearmesh cycles) %%
    if ~strcmp(inWheelDefectInfo.status,'Healthy') && strcmp(outWheelDefectInfo.status, 'Healthy')
        z2Repeat = inWheelInfo.z ;
    else
        z2Repeat = outWheelInfo.z ;
    end % of if
    
    gmsCyc = repmat(KeqOneGMCyc, z2Repeat, 1) ;
    gmsStruct.healthy = gmsCyc(1:simParam.gmsCoarseInterval:end) ;
    
    %% Build Keq struct for later calculations of the excitation force %%
    gmsStruct.KeqStruct.KeqHealthy = KeqHealthy ;
    gmsStruct.KeqStruct.gmCycDur = gmCycDur ;
    gmsStruct.KeqStruct.dCyc = simParam.dCycGMSFine ;
    gmsStruct.KeqStruct.z2Repeat = z2Repeat ;
    gmsStruct.KeqStruct.epsilon = inWheelInfo.epsilon ;
    
    %% Calculate and update gmsCyc for a damaged tooth %%
    % Clarification: The damage is seeded to the first tooth meshing,
    % meaning that the contribution of the first healthy tooth is removed
    % from gmsCyc and replaced by the stiffness of the damaged tooth.
    if ~strcmp(outWheelDefectInfo.status,'Healthy') || ~strcmp(inWheelDefectInfo.status,'Healthy')
        KeqDef = calcKeq(inWheelInfo, outWheelInfo, inWheelDefectInfo, outWheelDefectInfo) ;
        KeqDef = interp1(cycPtsCoarse, KeqDef, cycPtsFine, 'linear') ;
        gmsCyc(1:length(KeqDef)) = gmsCyc(1:length(KeqDef)) - KeqHealthy + KeqDef ; % update gmsCyc
        gmsStruct.KeqStruct.KeqFaulty = KeqDef ;
    end %% of if
    
    %% Generating the fields in the gms structure %%
    gmsStruct.gmsCyc = gmsCyc(1:simParam.gmsCoarseInterval:end) ;
    gmsStruct.dCyc = simParam.dCycGMSCoarse ;
    gmsStruct.cycVctr = (0:length(gmsStruct.gmsCyc)-1)' * simParam.dCycGMSCoarse ;
    gmsStruct.initContInd = outWheelInfo.initContInd ;
    
end % of function 'calcGearMeshStiff'