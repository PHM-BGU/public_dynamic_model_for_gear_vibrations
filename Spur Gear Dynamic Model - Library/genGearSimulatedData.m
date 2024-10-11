function [vibSigs, wheelsInfo, defectInfo, gmsStruct, ...
    EulerLagrangeComponents, numSol, testRigParam, simParam, ...
    operatingConds, initialConds] = genGearSimulatedData( ...
    module, Z, toothWidth, surfQuality, speedIn, loadOut, varargin)
%{
% Description:
% This core function drives the simulation model, generating synthetic
% vibration signals. It accepts essential parameters as input and
% systematically executes each step of the model's simulation process.
% =====
% Inputs:
% * module [mm] - the module of the transmission.
% * Z - a structure with the number of teeth on the wheels.
% * toothWidth [mm] - the width of the teeth on both wheels.
% * surfQuality - the precision grade of the teeth according to
%   DIN-3962 standard.
% * speedIn [rps] - the rotational speed on the driving shaft.
% * loadOut [N-m] - the torque applied to the output shaft.
% * The user can add optional inputs that will not get the default value
%   according to the dictionary below (buildDictionary) as follows:
%   get_patameters(..., 'Attribute', value)
% =====
% Outputs:
% * vibSigs - A structure with the simulated vibration signals and more
%   vital outputs for post analysis.
% * wheelsInfo - A structure with all the information about the wheels
%   (in a healthy status) and the manufacturing errors.
% * defectInfo - A structure with all the information about the faults.
% * gmsStruct - A structure with the gearmesh stiffness in cycle.
% * EulerLagrangeComponents - A structure holding the mass (M), damping (C)
%   and stiffness (K_cyc) matrices, as well as modal analysis.
% * numSol - The numerical solution, including displacement, velocity and
% acceleration matrices for each DOF.
% * [testRigParam, operatingConds, initialConds] - See getParam.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
    %}
    
    %% Generate structures with the model parameters %%
    [inWheelInfo, outWheelInfo, inWheelDefectInfo, outWheelDefectInfo, ...
        testRigParam, simParam, operatingConds, initialConds] = ...
        getParam(module, Z, toothWidth, surfQuality, speedIn, loadOut, varargin) ;
    
    %% Generate teeth profile & contact line properties %%
    [inWheelInfo.X, inWheelInfo.Y, inWheelInfo.R, inWheelInfo.toothModification.tipRelief.edgeInd] = ...
        genHealthyToothProfile(inWheelInfo.module, inWheelInfo.z, inWheelInfo.alpha, inWheelInfo.toothModification) ;
    [outWheelInfo.X, outWheelInfo.Y, outWheelInfo.R, outWheelInfo.toothModification.tipRelief.edgeInd] = ...
        genHealthyToothProfile(outWheelInfo.module, outWheelInfo.z, outWheelInfo.alpha, outWheelInfo.toothModification) ;
    [inWheelInfo, outWheelInfo] = calcContLineProp(inWheelInfo, outWheelInfo) ;
    
    %% Generate involute profile manufacturing errors %%
    inWheelInfo = genSurfQualityErrs(inWheelInfo) ;
    outWheelInfo = genSurfQualityErrs(outWheelInfo) ;
    [wheelsInfo.in, wheelsInfo.out] = deal(inWheelInfo, outWheelInfo) ;
    
    %% Generate tooth faults %%
    [outWheelDefectInfo, inWheelDefectInfo] = genDefectInfoStruct(...
        outWheelInfo, outWheelDefectInfo, inWheelInfo, inWheelDefectInfo) ;
    [defectInfo.in, defectInfo.out] = deal(inWheelDefectInfo, outWheelDefectInfo) ;
    
    %% Calculate gear mesh stiffness (gms) structure %%
    gmsStruct = calcGearMeshStiff(...
        inWheelInfo, outWheelInfo, simParam, ...
        inWheelDefectInfo, outWheelDefectInfo) ;
    
    %% Calculate Euler-Lagrange components %%
    EulerLagrangeComponents = calcEulerLagrangeComponents(...
        inWheelInfo, outWheelInfo, gmsStruct, testRigParam, ...
        simParam, operatingConds,inWheelDefectInfo, outWheelDefectInfo) ;
    
    %% Numerical solution %%
    [u, u_dot, u_ddot] = numericalSolution(EulerLagrangeComponents, initialConds, simParam) ;
    [numSol.disp, numSol.vel, numSol.acc] = deal(u, u_dot, u_ddot) ;
    
    %% Build a structure with the simulated signals %%
    % u = [xOut, yOut, xIn, yIn, thetaOut, thetaIn, thetaBrake, zOut, zIn, phiOut, phiIn, psiOut, psiIn]
    vibSigs.t = [0:length(u_ddot)-1]' / simParam.Fs ;
    vibSigs.hrz = (u_ddot(:, 1) + u_ddot(:, 3)) / simParam.g ; % from [m/s^2] to [g]
    vibSigs.vrt = (u_ddot(:, 2) + u_ddot(:, 4)) / simParam.g ; % from [m/s^2] to [g]
    vibSigs.axl = (u_ddot(:, 8) + u_ddot(:, 9)) / simParam.g ; % from [m/s^2] to [g]
    vibSigs.thetaDot.in = + u_dot(:, 6) / (2*pi) ; % from [r/s] to [Hz]
    vibSigs.thetaDot.out = - u_dot(:, 5) / (2*pi) ; % from [r/s] to [Hz]
    vibSigs.profErrMeshForceT = EulerLagrangeComponents.profErrMeshForcesT ;
    vibSigs.Fs = simParam.Fs ;
    vibSigs.operatingConds = operatingConds ;
    [vibSigs.units.time, vibSigs.units.sig] = deal('s', 'g') ;
    
end % of function 'genGearSimulatedData'