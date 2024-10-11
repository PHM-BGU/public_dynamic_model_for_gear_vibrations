function FexMtxT = calcForcesVctrT(inWheelInfo, outWheelInfo, cycMotorT, ...
    profErrMeshForceT, torqueBrakeT, torqueMotorT, inWheelDefectInfo, outWheelDefectInfo)
%{
% Description:
% This function calculates the excitation, or generalized external forces
% vector as a 2D matrix in time (ndof x lenT). The vector of coordinates:
% u = [xOut, yOut, xIn, yIn, thetaOut, thetaIn, thetaBrake, zOut, zIn, phiOut, phiIn, psiOut, psiIn]
% =====
% Inputs:
% * (in/out)WheelInfo - a structure with geomterical information of the tooth.
% * cycMotorT - [cyc], the cycle of the motor in time.
% * profErrMeshForceT - [N], a vector of the profile error induced forces,
%   obtained by multiplying the gms with the profile erros (displacement input).
% * torqueMotorT - [Nm], a vector of torsional torque in reaction to the
%   rotational speed of the motor.
% * torqueBrakeT - [Nm], a vector of torsinal load applied by the brake.
% * (in/out)WheelDefectInfo - a structure with all the information about the fault.
% =====
% Outputs:
% * FexMtxT - ndof x lenT matrix of the excitation force vector in time.
% =====
% Significant In-Function Variables:
% * (in/out)Rb - [m], base radius.
% * ndof - Number degrees of freedom ( = 13).
% * defectInfo.alphaCyc - [rad], a vector of the pressure angle in cycle
%   in case of a ThroughFaceFault.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
    %}
    
    %% Extract relevant parameters, mainly geometrical %%
    [alpha, beta] = deal(outWheelInfo.alpha, outWheelInfo.beta) ;
    [outR, inR] = deal(outWheelInfo.R, inWheelInfo.R)  ;
    [lenT, ndof] = deal(length(cycMotorT), 13) ;
    
    %% Calculate the variable Fex according to the health status %%
    switch outWheelDefectInfo.status
        case 'TroughFaceFault'
            FexMtxT = zeros(ndof, lenT) ;
            for ii = 1:lenT
                alphaT = outWheelDefectInfo.alphaT(ii) ;
                [geomZIndep, ~] = calcGeomCoeffMtx(alphaT, beta, outR, inR) ;
                FexMtxT(:, ii) = profErrMeshForceT .* geomZIndep ;
            end % of for ii
        otherwise
            [geomZIndep, geomZDep] = calcGeomCoeffMtx(alpha, beta, outR, inR) ;
            geomMtx = repmat(geomZIndep, 1, lenT) ;
            FexMtxT = geomMtx .* repmat(profErrMeshForceT', ndof, 1) ;
            cycT = rem(cycMotorT/(2*pi), outWheelInfo.z/inWheelInfo.z) ; % of the defected output wheel
            switch outWheelDefectInfo.status
                case 'ToothBreakage'
                    faultInds = find(cycT>=outWheelDefectInfo.defEntranceCycIn & cycT<=outWheelDefectInfo.defExitCycIn) ;
                    defWidth = outWheelDefectInfo.dimensions.z * 1e-3 ; % from [mm] to [m]
                    for ii = 1:length(faultInds)
                        zMean = -0.5*defWidth * ...
                            (outWheelDefectInfo.defExitCycIn - cycT(faultInds(ii))) / ...
                            (outWheelDefectInfo.defExitCycIn - outWheelDefectInfo.defEntranceCycIn) ;
                        geomVctr = geomZIndep + zMean*geomZDep ;
                        FexMtxT(:, faultInds(ii)) = profErrMeshForceT(faultInds(ii)) .* geomVctr ;
                    end % of for ii
                case 'PartialFaceFault'
                    faultInds = find(cycT>=outWheelDefectInfo.faultCycs.start & cycT<=outWheelDefectInfo.faultCycs.end) ;
                    zMean = -0.5*outWheelDefectInfo.dimensions.W * 1e-3 ; % from [mm] to [m]
                    geomVctr = geomZIndep + zMean*geomZDep ;
                    for ii = 1:length(faultInds)
                        FexMtxT(:, faultInds(ii)) = profErrMeshForceT(faultInds(ii)) .* geomVctr ;
                    end % of for ii
            end % of switch outWheelDefectInfo.status
    end % of switch outWheelDefectInfo.status
    
    FexMtxT(6, :) = FexMtxT(6, :) + torqueMotorT.' ; % adding the torque reaction from the motor
    FexMtxT(7, :) = FexMtxT(7, :) + torqueBrakeT.' ; % adding the external torque from the brake
    
end % of function 'calcForcesVctrT'
%%
function [geomZIndep, geomZDep] = calcGeomCoeffMtx(alpha, beta, outR, inR)
%{
% Description:
% This function calculates a matrix of geometrical coefficients of
% components in the stiffness matrix that are pressure-line-independent.
% This separation is meant to improve efficiency and running time.
% More information is detailed in the manual.
% =====
% Inputs:
% * alpha - [rad], pressure line (nominal or instantaneous).
% * beta - [rad], helix angle (zero for spur gears).
% * (in/out)R - [mm], structure with all radii.
% =====
% Outputs:
% * geomZ(In)Dep - an ndof x 1 vector of coefficients.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
%% Build the pressure-line-independent geometry vector %%
geomZIndep = [ ...
    - cos(beta)*sin(alpha), ...
    - cos(beta)*cos(alpha), ...
    + cos(beta)*sin(alpha), ...
    + cos(beta)*cos(alpha), ...
    - cos(beta)*outR.base * 1e-3, ... % from [mm] to [m]
    - cos(beta)*inR.base * 1e-3, ... % from [mm] to [m]
    0, ...
    - sin(beta), ...
    + sin(beta), ...
    0, ...
    0, ...
    + sin(beta)*outR.pitch * 1e-3, ... % from [mm] to [m]
    + sin(beta)*inR.pitch * 1e-3, ... % from [mm] to [m]
    ].' ;

%% Build the pressure-line-dependent geometry vector %%
geomZDep = [...
    zeros(1,9), ...
    + cos(alpha), ...
    - cos(alpha), ...
    - sin(alpha), ...
    + sin(alpha), ...
    ].' ;
end % of function 'calcGeomCoeffMtx'