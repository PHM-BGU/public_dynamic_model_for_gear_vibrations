function [inWheelInfo, outWheelInfo, inWheelDefectInfo, outWheelDefectInfo, ...
    testRigParam, simParam, operatingConds, initialConds] = ...
    getParam(module, Z, toothWidth, surfQuality, speedIn, loadOut, varargin)
%{
% Description:
% This function gets as an input the inspected transmission parameters,
% the desired test rig parameters, and the operational conditions, and
% returns structures with all the required information for the simulation.
% =====
% Inputs:
% * module [mm] - the module of the transmission.
% * Z - the number of teeth on the wheels:
%       # Z.in - number of teeth on the driving input wheel (pinion).
%       # Z.out - number of teeth on the driven output wheel (gear).
% * toothWidth [mm] - the width of the teeth on both wheels.
% * surfQuality - the precision grade of the teeth according to
%   DIN-3962 standard.
% * speedIn [rps] - the rotational speed on the driving shaft.
% * loadOut [Nm] - the torque applied to the output shaft.
% * The user can add optional inputs that will not get the default value
%   according to the dictionary below (buildDictionary) as follows:
%   getParam(..., 'Attribute', value)
% =====
% Outputs:
% * (in/out)WheelInfo - basic information about the wheel (in/out).
% * (in/out)WheelDefectInfo - fault parameters.
% * testRigParam - parameters of the test apparatus.
% * simParam - simulation parameters.
% * operatingConds - operational conditions (speed and load).
% * initialConds - initial conditions (displacement, velocity).
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
    %}
    
    %% Get operational conditions %%
    operatingConds.speed.rpsNom = speedIn ; % [Hz]
    operatingConds.speed.statisticErr = 0 ; % In default, no error
    operatingConds.load.loadNom = loadOut ; % [N-m]
    operatingConds.load.statisticErr = 0 ; % In default, no error
    
    %% Get initial conditions %%
    % u = [xOut, yOut, xIn, yIn, thetaOut, thetaIn, thetaBrake, zOut, zIn, phiOut, phiIn, psiOut, psiIn]
    ndof = 13 ; % ndof - Number Degrees of Freedom
    [initialConds.displacement, initialConds.velocity] = deal(zeros(ndof, 1)) ;
    initialConds.velocity(5:7) = speedIn*2*pi*[-Z.in/Z.out, 1, -Z.in/Z.out]' ;
    
    %% Get simulation parameters %%
    nominalFs = 1e5 ; % In default, the nominal sampling rate is 100kHz.
    simParam.numSol.trimPercent = 0.5 ;
    simParam.gmsCoarseInterval = 50 ; % Sample intervals for downsampling of the gear mesh stiffness function.
    simParam.numSol.gamma = 0.50 ; % A fixed value based on Newmark method.
    simParam.numSol.beta = 0.25 ; % A fixed value based on Newmark method.
    simParam.numSol.convergenceCriterion = 1e-5 ; % Convergence creterion for the numrical solution.
    simParam.numSol.maxIterationsNR = 10 ; % An early stopping callback of the maximum numer of iterations in the numerical solution.
    simParam.numSol.waitBar = false ; % A flag for displaying a waitbar during the numerical solution.
    simParam.g = 9.81 ; % [m/s^2]
    
    %% Fill test rig parameters %%
    testRigParam.eccentricity = 0 ; % [mm], Eccentricity level.
    testRigParam.zeta = 0.05 * ones(ndof, 1) ; % Damping ratio
    testRigParam.bearingStiff = 1.9e8 ; % [N/m]
    testRigParam.JBrake = 4.5e-3 ; % [kg-m^2]
    testRigParam.shaftDims.in.L = 200 ; % [mm]
    testRigParam.shaftDims.out.L = 600 ; % [mm]
    [testRigParam.shaftDims.in.R, testRigParam.shaftDims.out.R] = deal(10) ; % [mm]
    testRigParam.ndof = ndof ;
    
    %% Fill mutual info for both wheels %%
    wheelInfo.type = 'spur' ;
    wheelInfo.material = 'steel' ;
    wheelInfo.W = toothWidth ; % [mm]
    wheelInfo.module = module ; % [mm]
    wheelInfo.profile.initAngle = deg2rad(360) ; % [rad]
    wheelInfo.profile.NPts.profile = 10e3 ;
    wheelInfo.profile.NPts.invlt = 15e3 ;
    wheelInfo.profile.NPts.addendum = 5e3 ;
    wheelInfo.profile.NPts.fillet = 5e3 ;
    wheelInfo.alpha = deg2rad(20) ; % [rad], Pressure angle.
    wheelInfo.beta = deg2rad(0) ; % [rad], Base helix angle.
        wheelInfo.rotStiff = 1e8 ; % [N-m/rad], Rotational stiffness of the wheel.

    % Tooth profile modifications (tip relief & crowning)
    wheelInfo.toothModification.tipRelief.include = true ;
    wheelInfo.toothModification.tipRelief.shape = 'parabolic' ;
    wheelInfo.toothModification.tipRelief.xRelief = 1 ; % [mm], relieved tip loss (along x-axis).
    wheelInfo.toothModification.tipRelief.yRelief = 1 ; % [mm], relieved tip loss (along y-axis).
    wheelInfo.toothModification.crowning.include = false ;
    wheelInfo.toothModification.crowning.longAmount = 0.2 ; % [mm], longitudinal amount of crowningr(along y-axis).
    
    % Surface quality structure
    wheelInfo.surfQuality = surfQuality ;
    wheelInfo.surfQualityErrs.NPts = 1e2 ;
    wheelInfo.surfQualityErrs.statErrTyp = 0.1 ; % 10% of the nominal values
    
    % The coefficient matrix for the fillet foundation stiffness calculation
    coeffA = [-5.574 60.111 -50.952 -6.2042]*1e-5 ; %
    coeffB = [-1.9986 28.100 185.50 9.0889]*1e-3 ;
    coeffC = [-2.3015 -83.431 0.0538 -4.0964]*1e-4 ;
    coeffD = [4.77021 -9.9256 53.300 7.8297]*1e-3 ;
    coeffE = [0.0271 0.1624 0.2895 -0.1472] ;
    coeffF = [6.8045 0.9086 0.9236 0.6904] ;
    wheelInfo.coeffMtx4FilletStiff = [coeffA ; coeffB ; coeffC ; coeffD ; coeffE ; coeffF].' ;
    
    [inWheelInfo, outWheelInfo] = deal(struct) ;
    [inWheelDefectInfo.status, outWheelDefectInfo.status] = deal('Healthy') ; % In default, both wheels are healthy.
    
    %% Update variable argument inputs by parsing the optional inputs %%
    dict = buildDictionary ;
    
    if length(varargin) == 1 % in case the varargin are transfered across functions
        varargin = varargin{1} ;
    end % of if
    
    while ~isempty(varargin)
        key = varargin{1} ;
        variable = dict.(key) ;
        eval([variable ' = varargin{2};'])
        simParam.varargin.(key) = varargin{2} ;
        varargin(1:2) = [] ;
    end % of while
    
    %% Create info structure for each wheel %%
    [wheelInfo.E, wheelInfo.v, wheelInfo.G, wheelInfo.rho] = ...
        getMaterialProp(wheelInfo.material) ;
    
    for wheel = {'in', 'out'}
        name = wheel{1} ;
        infoStruct = wheelInfo ;
        infoStruct.wheel = name ;
        infoStruct.z = Z.(name) ; % Number of teeth.
        infoStruct = calcMassProp(infoStruct) ;
        infoStruct.shaftLen = testRigParam.shaftDims.(name).L ; % [mm]
        infoStruct.holeRadius = testRigParam.shaftDims.(name).R ; % [mm]
        
        try % in case of a desired profile errors vector
            infoStruct.surfQualityErrs = surfQualityErrsStruct.(name) ;
        end % of try
        
        if sum(size(testRigParam.zeta)) == 2 % uniform damping ratio
            testRigParam.zeta = testRigParam.zeta * ones(ndof, 1) ;
        elseif length(testRigParam.zeta) > ndof
            warning('only the first ndof=13 elements in zeta will be considered')
        end % of if
        
        testRigParam.shaftPolarStiff.(name) = calcShaftStiff(...
            testRigParam.shaftDims.(name).R, testRigParam.shaftDims.(name).L, wheelInfo.G) ;
        eval(sprintf('%sWheelInfo = infoStruct;', name)) ;
        clear infoStruct
    end % of for ii = 1:length(wheels)
    
    inWheelDefectInfo.statusOtherWheel = outWheelDefectInfo.status ;
    outWheelDefectInfo.statusOtherWheel = inWheelDefectInfo.status ;
    inWheelDefectInfo.zOtherWheel = outWheelInfo.z  ;
    outWheelDefectInfo.zOtherWheel = inWheelInfo.z ;
    [inWheelDefectInfo.zInWheel, outWheelDefectInfo.zInWheel] = deal(Z.in) ;
    [inWheelDefectInfo.NPtsY4Shear, outWheelDefectInfo.NPtsY4Shear] = deal(1e4) ;
    [inWheelDefectInfo.dimensions.circNPts, outWheelDefectInfo.dimensions.circNPts] = deal(1e4) ;
    
    %% Get simulation paramters %%
    lcmTeeth = lcm(Z.in, Z.out) ;
    simParam.dCycGMSFine = ...
        1 / ( ...
        simParam.gmsCoarseInterval * lcmTeeth * ...
        ceil (max(nominalFs, 1e6) / simParam.gmsCoarseInterval / lcmTeeth)) ;
    simParam.dCycGMSCoarse = simParam.dCycGMSFine * simParam.gmsCoarseInterval ;
    simParam.Fs = ceil(nominalFs/(lcmTeeth * speedIn)+6) * lcmTeeth * speedIn ;
    simParam.NumInCycs = 2 * Z.out ; % of the input driving wheel
    simParam.numSol.sigFirstInd = (simParam.Fs / speedIn) * ...
        floor(simParam.numSol.trimPercent * simParam.NumInCycs) + 1 ;
    
    %% Allow variable rotational speed %%
    try
        operatingConds.speedInSigT = speedInSigT ; % here 'speedInSigT' is a varargin.
    catch
        T = simParam.NumInCycs / speedIn ;
        dt = 1/simParam.Fs ;
        t = (0:dt:(T-dt))' ;
        speedInSigT = uniformRand(speedIn, operatingConds.speed.statisticErr, size(t)) ;
        operatingConds.speedInSigT = speedInSigT ;
    end % of try
    
end % of function 'getParam'

%%
function dict = buildDictionary()
%% Attributes of the wheels %%
dict.Alpha = 'wheelInfo.alpha' ; % [rad], Pressure angle.
dict.RotStiff = 'wheelInfo.rotStiff' ; % Rotational stiffness of the wheel [N-m/rad]
dict.Material = 'wheelInfo.material' ; % The wheel's metal
dict.InitialProfileAngle = 'wheelInfo.profile.initAngle' ; % [rad], Initial profile angle.
dict.NptProfile = 'wheelInfo.profile.NPts.profile' ;
dict.NptInvlt = 'wheelInfo.profile.NPts.invlt' ;
dict.NptAddendum = 'wheelInfo.profile.NPts.addendum' ;
dict.NptFillet = 'wheelInfo.profile.NPts.fillet' ;
dict.MassProps = 'wheelInfo.massProps' ;
dict.TipReliefStruct = 'wheelInfo.tipRelief' ; % a structure with all necessary fields for generating tip relief
%% Attributes of the profile errors %%
dict.SurfQualityErrs = 'surfQualityErrsStruct' ; %
dict.NPts4SurfQuality = 'wheelInfo.surfQualityErrs.NPts' ;  % Number of points for low resolution (coarse) in the profile error function
dict.StatErr4SurfQ = 'wheelInfo.surfQualityErrs.statErrTyp' ; % Typical percentage of the nominal value in the profile error parameters for scattering
dict.CoeffMtx4FilletStiff = 'wheelInfo.coeffMtx' ; % Coefficient matrix for fillet foundation matrix.
%% Attributes of the fault %%
dict.InFaultProp = 'inWheelDefectInfo' ; % Properties of the fault of the input driving wheel
dict.OutFaultProp = 'outWheelDefectInfo' ; % Properties of the fault of the output driven wheel
dict.NPts4Shear = 'outWheelDefectInfo' ; % Property of certain faults in the calculation of the shear potential energy.
%% Attributes of the test apparatus %%
dict.Eccentricity = 'testRigParam.eccentricity' ; % The level of eccentricity [m]
dict.Zeta = 'testRigParam.zeta' ; % damping ratio
dict.BearingStiff = 'testRigParam.bearingsStiff' ; % [N/m], Typical stiffness of the support bearings.
dict.BrakeInertia = 'testRigParam.JBrake' ; % [kg-m^2], Mass moment of inertia of the loading source ('brake').
dict.InShaftLen = 'testRigParam.shaftDims.in.L' ; % [mm], Input shaft length between the motor and the pinion.
dict.OutShaftLen = 'testRigParam.shaftDims.out.L' ; % [mm], Output shaft length between the brake and the gear.
dict.InShaftRadius = 'testRigParam.shaftDims.in.R' ; % [mm], Radius of the input shaft (aka the input wheel's hole).
dict.OutShaftRadius = 'testRigParam.shaftDims.out.R' ; % [mm], Radius of the output shaft (aka the output wheel's hole).
%% Attributes of the resolution and numerical solution %%
dict.Intervals4GMS = 'simParam.gmsCoarseInterval' ; % Sample intervals for downsampling of the gear mesh stiffness function.
dict.SamplingRate = 'nominalFs' ; % [Hz or S/s], The nominal sampling rate of the signal.
dict.ConvergenceCriterion = 'simParam.numSol.convergenceCriterion' ; % Convergence creterion for the numrical solution.
dict.MaxIterations = 'simParam.numSol.maxIterationsNR' ; % An early stopping callback of the maximum numer of iterations in the numerical solution.
dict.NumericalSolMethod = 'simParam.numSolMethod' ; % The method of the numerical solution.
dict.SigTrimPercent = 'simParam.numSol.trimPercent' ; % The percentage of data to trim from the signal (due to numerical convergence).
dict.WaitBar = 'simParam.numSol.waitBar' ; % A flag to showing waitbar during the numerical solution (0-no, else-yes).
%% Attributes of the operational conditions %%
dict.SpeedStatError = 'operatingConds.speed.statisticErr' ; % [rps], Statistical error of the speed (U(-a,a) where a is the parameter).
dict.LoadStatError = 'operatingConds.load.statisticErr' ; % [N-m], Statistical error of the load (U(-a,a) where a is the parameter).
dict.InitialConditionDisp = 'initialConds.displacement' ; % Vector of the initial conditions of the displacement.
dict.InitialConditionVel = 'initialConds.velocity' ; % Vector of the initial conditions of the velocity.
dict.SpeedInSigT = 'speedInSigT' ; % A desired shaft rotational speed profile.
end % of function 'buildDictionary'