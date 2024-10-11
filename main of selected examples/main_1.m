%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Dynamic Model of Spur Gears  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%            main_1 - Healthy             %%%%%%%%%%%%%
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
[module, Z.in, Z.out, toothWidth] = deal(3, 17, 38, 15) ; % (3, 18, 35, 30) ;
[speedIn, loadOut, surfQuality] = deal(40, 10, 7) ; % (45, 15, 0) ;
[FsDownSamp, chan] = deal(25e3, 'vrt') ; 

%% Generate simulation %%
tic
[vibSigs, wheelsInfo, defectInfo, gmsStruct, EulerLagrangeComponents, ...
    numSol, testRigParam, simParam, operatingConds, initialConds] = ...
    genGearSimulatedData( ...
    module, Z, toothWidth, surfQuality, speedIn, loadOut, ...
    'SamplingRate', 1e5, 'WaitBar', true) ;
toc

%% Signal processing %%
saStruct = buildSAStruct(vibSigs.(chan), speedIn*Z.in/Z.out, Z.out, simParam.Fs, FsDownSamp) ;

%% Display simulation results %%
dispSimResults(vibSigs, wheelsInfo, defectInfo, gmsStruct, 'VibSigChannel', chan)