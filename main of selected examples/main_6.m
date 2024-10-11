%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Dynamic Model of Spur Gears  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%      main_6 - Tooth Destruction        %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% © PHM-BGU Laboratory, Ben-Gurion University of the Negev, %%%%
%%%%%%%%%%%%%%%%   Be'er Sheva, Israel. 2023     %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; clc ; close all ;
path2fold = 'C:\Users\Lior\OneDrive - post.bgu.ac.il\Project - OL - 2023\Model API' ;
addpath(fullfile(path2fold, 'main of selected examples'))
addpath(fullfile(path2fold, 'Spur Gear Dynamic Model - Library'))
addpath(fullfile(path2fold, 'Gear Signal Processing - Library'))

%% Define basic parameters of the desired apparatus %%
[module, Z.in, Z.out, toothWidth] = deal(3, 18, 35, 30) ; % deal(3, 17, 38, 15)
[speedIn, loadOut, surfQuality] = deal(45, 15, 8) ;
[FsDownSamp, chan] = deal(50e3, 'vrt') ;

%% Define parameters of the desired fault %%
defectInfo.status = 'ToothDestruction' ;
defectInfo.dimensions.method = 'byDefEndCoeff' ;
defectInfo.dimensions.defStartCoeff = 3/14 ;
defectInfo.dimensions.yOffset = 0 ; % [mm]
defectInfo.dimensions.noiseLevel = 0 ; % [mm]

switch defectInfo.dimensions.method
    case 'byDefEndCoeff'
        defectInfo.dimensions.defEndCoeff = 13/14 ;
    case 'byDefSlope'
        defectInfo.dimensions.defSlope = deg2rad(22) ; % from [deg] to [rad]
    otherwise
        error('The chosen method for fault dimensions does not exist')
end % of switch-case

%% Generate simulation %%
tic
[vibSigs, wheelsInfo, defectInfo, gmsStruct, EulerLagrangeComponents, ...
    numSol, testRigParam, simParam, operatingConds, initialConds] = ...
    genGearSimulatedData( ...
    module, Z, toothWidth, surfQuality, speedIn, loadOut, ...
    'OutFaultProp', defectInfo, ...
    'SamplingRate', 1e5, 'WaitBar', true) ;
toc

%% Signal processing %%
saStruct = buildSAStruct(vibSigs.(chan), speedIn*Z.in/Z.out, Z.out, simParam.Fs, FsDownSamp) ;

%% Display simulation results %%
dispSimResults(vibSigs, wheelsInfo, defectInfo, gmsStruct, 'VibSigChannel', chan)