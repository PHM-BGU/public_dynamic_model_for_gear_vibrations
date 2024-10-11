%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Dynamic Model of Spur Gears  %%%%%%%%%%%%%%%%%%%
%%%   main_7 - Multiple Faults With the Same Surface Quality    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% © PHM-BGU Laboratory, Ben-Gurion University of the Negev, %%%%
%%%%%%%%%%%%%%%%   Be'er Sheva, Israel. 2023     %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; clc ; close all ;
path2fold = 'C:\Users\Lior\OneDrive - post.bgu.ac.il\Project - OL - 2023\Model API' ;
addpath(fullfile(path2fold, 'main of selected examples'))
addpath(fullfile(path2fold, 'Spur Gear Dynamic Model - Library'))
addpath(fullfile(path2fold, 'Gear Signal Processing - Library'))
path2save = fullfile(path2fold, 'SimulatedData') ;

%% Define basic parameters of the desired apparatus %%
[module, Z.in, Z.out, toothWidth] = deal(3, 18, 35, 30) ; % deal(3, 17, 38, 15)
[speedIn, loadOut, surfQuality] = deal(45, 15, 8) ;
[FsDownSamp, chan] = deal(50e3, 'vrt') ;

%% Generate a healthy simulation & signal processing %%
tic
[vibSigs, wheelsInfo] = genGearSimulatedData( ...
    module, Z, toothWidth, surfQuality, speedIn, loadOut, ...
    'SamplingRate', 1e5, 'WaitBar', false) ;
toc

figure(1) ; hold on ;
plot(vibSigs.t, vibSigs.(chan), 'DisplayName', 'H') ;
saStruct = buildSAStruct(vibSigs.(chan), speedIn*Z.in/Z.out, Z.out, vibSigs.Fs, FsDownSamp) ;
save([path2save, '\healthy'], 'saStruct') ;
vibSigsH = vibSigs ;

%% Extract the surface quaility erros from the healthy simulation %%
SurfQualityErrs.in = wheelsInfo.in.surfQualityErrs ;
SurfQualityErrs.out = wheelsInfo.out.surfQualityErrs ;
save([path2save, '\SurfQualityErrs'], 'SurfQualityErrs') ;
clear vibSigs wheelsInfo saStruct

%% Define parameters of the desired fault %%
defectInfoStatus = 'ToothDestruction' ;

switch defectInfoStatus
    case 'ThroughFaceFault'
        defectInfoDimMethod = 'byCircCoordinates' ;
        xCens = [57.12, 57.69, 57.37, 57.49] ;
        yCens = [3.59, 3.71, 4, 3.66] ;
        rCircs = [1.7, 2.1, 2.4, 2.4] ;
        numFaults = length(xCens) ;
    case 'ToothDestruction'
        defectInfoDimMethod = 'byDefEndCoeff' ;
        defStartCoeffs = [+05.0, +03.5, +00.2, -00.6, -00.8, -01.0] / 14 ;
        defEndCoeffs =   [+12.0, +12.2, +13.0, +13.2, +13.2, +13.8] / 14 ;
        yOffsets =       [000.0, 000.0, 000.0, 000.0, 000.0, 000.0] ;
        numFaults = length(defStartCoeffs) ;
    otherwise
        error('The chosen fault does not exist')
end % of switch-case

%% Pre-Allocation %%
vibSigsArr = cell(numFaults+1, 1) ;
vibSigsArr{1} = vibSigsH ;

%% Generate faulty simulation iteratively %%
for ii = 1 : numFaults
    defectInfo.status = defectInfoStatus ;
    switch defectInfoStatus
        case 'ThroughFaceFault'
            defectInfo.dimensions.method = defectInfoDimMethod ;
            [defectInfo.dimensions.circCen.x, ...
                defectInfo.dimensions.circCen.y, ...
                defectInfo.dimensions.circR] = deal(xCens(ii), yCens(ii), rCircs(ii)) ;
        case 'ToothDestruction'
            defectInfo.dimensions.method = defectInfoDimMethod ;
            defectInfo.dimensions.defEndCoeff = defEndCoeffs(ii) ;
            defectInfo.dimensions.defStartCoeff = defStartCoeffs(ii) ;
            defectInfo.dimensions.yOffset = yOffsets(ii) ;
    end % of switch-case
    
    tic
    vibSigs = genGearSimulatedData( ...
        module, Z, toothWidth, surfQuality, speedIn, loadOut, ...
        'OutFaultProp', defectInfo, ...
        'SurfQualityErrs', SurfQualityErrs, ...
        'SamplingRate', 1e5, 'WaitBar', false) ;
    toc
    
    h = plot(vibSigs.t, vibSigs.(chan), 'DisplayName', [defectInfoStatus, num2str(ii)]) ;
    uistack(h, 'bottom') ;
    vibSigsArr{1+ii} = vibSigs ;
    saStruct = buildSAStruct(vibSigs.(chan), speedIn*Z.in/Z.out, Z.out, vibSigs.Fs, FsDownSamp) ;
    save([path2save,'\', defectInfo.status, num2str(ii)], 'saStruct') ;
    clear vibSigs defectInfo saStruct
    
end % of for ii

%% Figure settings %%
set(gca,'FontName','Cambria','FontSize',10)
xlabel('Time [s]') ; ylabel('Amplitude [g]')
legend('Location', 'northwest')

%% Save results %%
save(fullfile(path2save, 'vibSigs'), 'vibSigsArr') ;