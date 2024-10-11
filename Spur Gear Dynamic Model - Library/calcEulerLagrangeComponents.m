function EulerLagrangeComponents = calcEulerLagrangeComponents(...
    inWheelInfo, outWheelInfo, gmsStruct, testRigParam, ...
    simParam, operatingConds, inWheelDefectInfo, outWheelDefectInfo)
%{
% Description:
% This function calculates the Mass (M), Damping (C), and Stiffness (K)
% matrices, and the external non-conservative force vector (Fex) according
% to Euler-Lagrange model for the equations of motion.
% Notice that both Fex and K are time-variable. Hence, the
% external force is a 2D mtx in time while K is a 3D array in cycle.
% In addition, the natural frequencies and the natural modes are
% calculated and documented.
% =====
% Inputs:
% * [inWheelInfo, outWheelInfo] - basic information about the wheels.
% * gmsStruct - a structure with the gms in cycle.
% * testRigParam - s structure with structural properties of the test rig:
%   # testRigParam.JBrake - [kg-m^2], Polar moment on intertia of the brake.
%   # testRigParam.zeta - a vector of modal damping ratios.
%   # testRigParam.eccentricity - [mm], eccentricity level.
%   # testRigParam.shaftPolarStiff.(in/out) - [Nm/rad], Polar stiffness of
%     torsional shafts (in and out).
%   # testRigParam.bearingStiff - [N/m], The stiffness of the springs
%     representing the support bearings (typical).
% * simParam - a structure with parameters related to the simulation:
%   # simParam.NumInCycs - duration of the simulated signal, specified in 
%     terms of the number of cycles completed by the driving wheel.
%   # simParam.Fs - [S/s], sampling rate.
%   # simParam.dCycGMSCoarse - [cyc], resolution of the gms in cycle.
%   # simParam.numSol.sigFirstInd - the first index in the simulated signal
%     after trimming, used in this function for documentation purposes.
% * operatingConds - A structure with nominal values and statistic errors
%   of the operating conditions, including input speed and output load.
% * in/outWheelDefectInfo - Structures with information about tooth faults.
% =====
% Outputs:
% * EulerLagrangeComponents - a structure with the follwing fields:
%   # M - generalized diagonal mass mtx [ndof x ndof]
%   # KCycList - variable stiffness mtx in cycle [ndof x ndof x cycle].
%   # C - modal damping mtx (assuming uniform damping ratio).
%   # modalMtx - natural vibration modes (eigenvectors).
%   # naturalFreqs - [Hz] - natural frequencies vector.
%   # FexVctrT - excitation forces vector in time [ndof x t] .
%   # alpha.time/cyc - [rad], documentation of the variable pressure angle.
%   # profErrMeshForcesT - [N], the mesh forces induced by profile errors.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
    %}
    
    %% Building the motor's cycle vector (in time) %%
    rpsNom = operatingConds.speed.rpsNom ;
    T = simParam.NumInCycs / rpsNom ;
    dt = 1 / simParam.Fs ;
    timeVctr = (0:dt:(T-dt))' ;
    rpsSig = uniformRand(rpsNom, operatingConds.speed.statisticErr, size(timeVctr)) ;
    cycMotor = 2*pi*cumsum(rpsSig*dt) ;
    cycMotor = cycMotor - cycMotor(1) ;
    tr = outWheelInfo.z / inWheelInfo.z ;
    gmsCyc = gmsStruct.gmsCyc * 1e3 ; % from [N/mm] to [N/m]
    
    %% Calculate the generalized mass mtx M %%
    M = calcGeneralizedMassMtx(inWheelInfo, outWheelInfo, testRigParam.JBrake) ;
    
    %% Calculate the stiffness mtx K (3D array in cycle) %%
    outWheelDefectInfo.alphaCyc = calcAlphaCyc(outWheelInfo.alpha, ...
        gmsStruct.cycVctr, tr, outWheelDefectInfo) ;
    KCycList = calcVariableStiffMtx(gmsCyc, gmsStruct.cycVctr, outWheelInfo, ...
        inWheelInfo, testRigParam, outWheelDefectInfo, inWheelDefectInfo) ;
    
    %% Solve eigenvalue problem and calculate the nautral frequencies %%
    avgK = sum(KCycList, 3) / length(gmsCyc)  ;
    [modalMtx, spectralMtx] = eig(avgK, M) ; % eigenvalues
    omega = sqrt(diag(spectralMtx)) ;
    
    %% Calculate the damping mtx C %%
    zeta = testRigParam.zeta ;
    ndof = length(omega) ;
    C = zeros(ndof) ;
    for n = 1:ndof
        modalMass = modalMtx(:,n)' * M * modalMtx(:,n) ;
        modalDamping = 2 * zeta(n) * omega(n) / modalMass ;
        C = C + M*(modalDamping*modalMtx(:,n)*(modalMtx(:,n))')*M ;
    end % of for n
    
    %% Calculate mesh forces induced by surface quality errors %%
    [inSurfQMeshForceT, outSurfQMeshForceT] = deal(zeros(size(cycMotor))) ;
    if inWheelInfo.surfQuality
        inSurfQMeshForceT = calcProfErrMeshForcesT( ...
            inWheelInfo.surfQualityErrs.errsMtx, gmsStruct.KeqStruct, 'spline', dt, rpsSig) ;
    end % of if
    %
    if outWheelInfo.surfQuality
        outSurfQMeshForceT = calcProfErrMeshForcesT( ...
            outWheelInfo.surfQualityErrs.errsMtx, gmsStruct.KeqStruct, 'spline', dt, rpsSig) ;
    end
    surfQMeshForceT = inSurfQMeshForceT + outSurfQMeshForceT ;
    
    %% Calculate mesh forces induced by eccentricity errors %%
    eccentricityMeshForceT = calcEccentricityMeshForceT( ...
        testRigParam.eccentricity, gmsStruct, tr, outWheelInfo.alpha, dt, rpsSig) ;
    
    %% Calculate mesh forces induced by fault displacement %%
    faultDispMeshForceT = calcFaultDispMeshForceT( ...
        outWheelDefectInfo, gmsStruct, dt, rpsSig, tr) ;
    
    %% Calculate the excitation force vector in time Fex %%
    profErrMeshForcesT = ...
        surfQMeshForceT + eccentricityMeshForceT + faultDispMeshForceT ;
    torqueMotorT = testRigParam.shaftPolarStiff.in * cycMotor ;
    torqueBrakeT = uniformRand(operatingConds.load.loadNom, operatingConds.load.statisticErr, size(timeVctr)) ;
    outWheelDefectInfo.alphaT = convertCyc2Time(simParam.dCycGMSCoarse, outWheelDefectInfo.alphaCyc, dt, rpsSig) ;
    FexVctrT = calcForcesVctrT(inWheelInfo, outWheelInfo, cycMotor, ...
        profErrMeshForcesT, torqueBrakeT, torqueMotorT, inWheelDefectInfo, outWheelDefectInfo) ;
    
    %% Organize all Euler-Lagrange components in a structure %%
    [EulerLagrangeComponents.M, EulerLagrangeComponents.C, ...
        EulerLagrangeComponents.KCycList, EulerLagrangeComponents.FexVctrT] = ...
        deal(M, C, KCycList, FexVctrT) ;
    EulerLagrangeComponents.naturalFreqs = omega / (2*pi) ; % from [r/s] to [Hz]
    EulerLagrangeComponents.modalMtx = modalMtx ;
    EulerLagrangeComponents.rpsSig = rpsSig(simParam.numSol.sigFirstInd:end) ;
    EulerLagrangeComponents.alpha.cyc = outWheelDefectInfo.alphaCyc ;
    EulerLagrangeComponents.alpha.time = outWheelDefectInfo.alphaT(simParam.numSol.sigFirstInd:end) ;
    EulerLagrangeComponents.profErrMeshForcesT = profErrMeshForcesT(simParam.numSol.sigFirstInd:end) ;
    
end % of function 'calcEulerLagrangeComponents'
%%
function M = calcGeneralizedMassMtx(inWheelInfo, outWheelInfo, brakeInertia)
%{
% Description:
% This function calculates the constant diagonal mass mtx, according to
% the following vector of generalized coordinates:
% u = [xOut, yOut, xIn, yIn, thetaOut, thetaIn, thetaBrake, zOut, zIn, phiOut, phiIn, psiOut, psiIn]
% =====
% Inputs:
% * wheelInfo - in/out wheel mass properties:
    # wheelInfo.m - [kg], mass.
    # wheelInfo.rotI - [kg-m^2] - rotational mass moment of inertia.
    # wheelInfo.JRot - [kg-m^2] - polar mass moment of inertia.
% * brakeInertia - [kg-m^2] - polar mass moment of inertia of the brake.
% =====
% Outputs:
% * M - The diagonal generalized mass mtx.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
M = diag([outWheelInfo.m, outWheelInfo.m, inWheelInfo.m, ...
    inWheelInfo.m, outWheelInfo.polarJ, inWheelInfo.polarJ, ...
    brakeInertia, outWheelInfo.m, inWheelInfo.m, outWheelInfo.rotI, ...
    inWheelInfo.rotI, outWheelInfo.rotI, inWheelInfo.rotI]) ;
end % of function 'calcGeneralizedMassMtx'
%%