function faultDispMeshForceT = calcFaultDispMeshForceT(defectInfo, gmsStruct, dt, rpsSig, tr)
%{
% Description:
% This function calculates the mesh forces induced by fault displacement.
% The displacement results from local or distributed tooth faults causing
% deviations in the involute tooth profile. The calculated displacement
% takes both wheels (in and out) into account.
% Meshing force is computed as the product of fault displacement and
% gear mesh stiffness, both over time.
% Note that the 'ThroughFaceFault' displacement requires direct
% time-dependent calculations, preventing the use of the 'calcMeshForcesT'
% function in its intended form. However, for faults like 'ToothDestruction'
% that don't necessitate time-dependent calculations, this function calls
% the parent function 'calcMeshForcesT'.
% =====
% Inputs:
% * defectInfo - provides fundamental information about the wheels. For
%   more detailed specifics, please refer to the description in 'genXXXProfile,'
%   where XXX represents the inspected fault.
% * gmsStruct - gearmesh stiffness structure. In the case of faults
%   requiring time-dependent calculations, the gearmesh stiffness
%   (gmsStruct.gmsCyc) is converted from cycles to time. Otherwise, a
%   designated structure (gmsStruct.KeqStruct) for meshing force calculations
%   is utilized. For more information, please refer to 'calcMeshForcesT.'
% * [rpsSig, dt], - input speed signal [Hz] and the time resolution [s].
% * tr - Transmission ratio, defined as the ratio between the number of
%   teeth, as follows: zOut/zIn.
% =====
% Outputs:
% * faultDispMeshForceT - [N], The mesh forces vector in time.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
switch defectInfo.status
    case 'ThroughFaceFault'
        cycMotor = 2*pi*cumsum(rpsSig*dt) ;
        cycMotor = cycMotor - cycMotor(1) ;
        [h, defStartCyc, defEndCyc] = deal( ...
            defectInfo.h, defectInfo.faultCycs.start, defectInfo.faultCycs.end) ;
        faultDisplacement = zeros(size(cycMotor)) ;
        for ii = 1:length(cycMotor)
            currentCyc = rem(cycMotor(ii)/(2*pi), tr) ;
            if currentCyc>defStartCyc && currentCyc<defEndCyc
                realtiveCyc = (currentCyc-defStartCyc) / (defEndCyc-defStartCyc) ;
                faultDisplacement(ii) = -h * (1-cos(2*pi*realtiveCyc)) / 2 ;
            end % of if
        end % of for ii
        gmsT = convertCyc2Time(gmsStruct.dCyc, gmsStruct.gmsCyc, dt, rpsSig) ;
        faultDispMeshForceT = faultDisplacement .* gmsT ;
    case 'ToothDestruction'
        faultDispMeshForceT = calcProfErrMeshForcesT(defectInfo.faultDisplacement, ...
            gmsStruct.KeqStruct, 'PCHIP', dt, rpsSig) ;
    otherwise
        faultDispMeshForceT = 0 ;
end % of switch defectInfo.status
end % of function 'calcFaultDispMeshForceT'