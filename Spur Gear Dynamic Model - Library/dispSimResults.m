function dispSimResults(vibSigs, wheelsInfo, defectInfo, gmsStruct, varargin)
% {
% Description:
% This function displayes a visual panel with the main simulation results.
% The panel is a 3x5 figure, divided as follows:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    |        1          |         2         |         3            |   4-5    |  %
%    ---------------------------------------------------------------------------  %
%  1 | Tooth in YX (In ) | Tooth in ZX (In ) | Profile Errors (In ) | SA (In ) |  %
%    ---------------------------------------------------------------------------  %
%  2 | Tooth in YX (Out) | Tooth in ZX (Out) | Profile Errors (Out) | SA (Out) |  %
%    ---------------------------------------------------------------------------  %
%  3 |                          Gearmesh Stiffness                             |  %
%    ---------------------------------------------------------------------------  %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====
% Inputs:
% * vibSigs - a structure with the vibration signals. The default channel
% is set to the vertical direction (vrt), but the user is allowed to modify
% it in the varargin.
% * wheelsInfo - a structure with all the wheels information (in/out).
% * defectInfo - a structure with all the fault related information (in/out).
% * gmsStruct - a structure with the gearmesh stiffness in cycle.
% * varargin - variable arguments in according to the dictionary defined at the bottom of this function.
% =====
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
% }

%% Default Parameters Settings %%
titleFlag = false ;
colorSurfQErrs = [0.6 0.4 0.08] ;
colorSA = [0.33, 0.35, 0.65] ;
colorGMSHealthy = [0.75 0.75 0.75] ;
vibSigChan = 'vrt' ;

%% Load Variable Argument Inputs %%
dict = buildDictionary ;
while ~isempty(varargin)
    key = varargin{1} ;
    variable = dict.(key) ;
    eval([variable ' = varargin{2};'])
    varargin(1:2) = [] ;
end % of while

figure
%% Tooth Profiles in YX Plane %%
inToothYX = subplot(3, 5, 1) ;
plotToothProfile(wheelsInfo.in, defectInfo.in) ;
title('Tooth Profile (In)') ; plotSettings ;
%
outToothYX = subplot(3, 5, 6) ;
plotToothProfile(wheelsInfo.out, defectInfo.out) ;
title('Tooth Profile (Out)') ; plotSettings ;

%% Tooth Faces in ZX Plane %%
inToothZX = subplot(3, 5, 2) ;
plotToothFace(wheelsInfo.in, defectInfo.in) ;
title('Tooth Face (In)') ; plotSettings ;
%
outToothZX = subplot(3, 5, 7) ;
plotToothFace(wheelsInfo.out, defectInfo.out) ;
title('Tooth Face (Out)') ; plotSettings ;

linkaxes([inToothYX inToothZX],'y')
linkaxes([outToothYX outToothZX],'y')

%% Manufacturing Profile Errors %%
inProfErrs = subplot(3, 5, 3) ;
if wheelsInfo.in.surfQuality
    surfQualErr2plot = wheelsInfo.in.surfQualityErrs.errsMtx(:,1) * 1e3 ;
    plot(surfQualErr2plot, 'color', colorSurfQErrs)
    set(gca, 'XTick', []);
    xticks([1, length(surfQualErr2plot)]);
    xticklabels({'initCont', 'Tip'});
else
    plot([0 1], [0 0], 'color', colorSurfQErrs)
end % of if
xlabel('First Tooth') ;  ylabel('Errors [\mum]') ;
title('Manufacturing Errors (In)') ; plotSettings ;
%
outProfErrs = subplot(3, 5, 8) ;
if wheelsInfo.out.surfQuality
    surfQualErr2plot = flip(wheelsInfo.out.surfQualityErrs.errsMtx(:,1)) * 1e3 ;
    plot(surfQualErr2plot, 'color', colorSurfQErrs)
    set(gca, 'XTick', []);
    xticks([1, length(surfQualErr2plot)]);
    xticklabels({'initCont', 'Tip'});
else
    plot([0 1], [0 0], 'color', colorSurfQErrs)
end % of if
xlabel('First Tooth') ;  ylabel('Errors [\mum]') ;
title('Manufacturing Errors (Out)') ; plotSettings ;

%% Synchronous Average (SA) %%
vibSig = vibSigs.(vibSigChan) ;
%
saInSig = mean(reshape(vibSig, [length(vibSig)/wheelsInfo.out.z, wheelsInfo.out.z]), 2) ;
cycIn = [0:length(saInSig)-1]'/length(saInSig) ;
saIn = subplot(3, 5, 4:5) ;
plot(cycIn, saInSig, 'color', colorSA)
xlabel('Cycle (in)') ;  ylabel('Amplitude [g]') ;
title('SA - In Shaft') ; plotSettings ;
%
saOutSig = mean(reshape(vibSig, [length(vibSig)/wheelsInfo.in.z, wheelsInfo.in.z]), 2) ;
cycOut = [0:length(saOutSig)-1]'/length(saOutSig) ;
saOut = subplot(3, 5, 9:10) ;
plot(cycOut, saOutSig, 'color', colorSA)
xlabel('Cycle (out)') ;  ylabel('Amplitude [g]') ;
title('SA - Out Shaft') ; plotSettings ;

%% Gearmesh Stiffness (GMS) in Cycle (In) %%
gmsInCycle = subplot(3, 5, 11:15) ;
if ~strcmp(defectInfo.out.status, 'Healthy')
    plot(gmsStruct.cycVctr, gmsStruct.healthy, '--', 'color', colorGMSHealthy)
    hold on ;
end % of if
plot(gmsStruct.cycVctr, gmsStruct.gmsCyc, 'color', [0 0.5 0.4])
xlim([0, gmsStruct.cycVctr(end)])
xlabel('Cycle (in)') ;  ylabel('gms [N/mm]') ;
title('Gear Mesh Stiffness') ; plotSettings ;

%% Set Title %%
try if titleFlag
        titleHandle = sgtitle({...
            ['Simulation Results of Gear Pair:  (zI_n, zO_uT, module)=', ...
            num2str(zIn),',',num2str(zOut),',',num2str(wheelsInfo.in.module),...
            ';   input speed=',num2str(inSpeed), ...
            'rps ;   output load=',num2str(outLoad), ...
            'Nm ;   surface quality= DIN', num2str(wheelsInfo.in.surfQuality)] ; ...
            ['Output Gear Status:  ', defectInfo.out.status]}) ;
        set(titleHandle, 'FontName', 'Cabmria', 'FontSize', 8, 'FontWeight', 'bold') ;
    end % of if
end % of try

end % of function 'dispSimResults'

%%
function plotSettings()
set(gca, 'FontName', 'Calibri','FontSize', 8) ;
box on ;
end % of function 'plotSettings'

%%
function dict = buildDictionary()
dict.VibSigChannel = 'vibSigChan' ; % desired channel to plot
dict.Title = 'titleFlag' ; % title flag
dict.ColorSurfQErrs = 'colorSurfQErrs' ;
dict.ColorSA = 'colorSA' ;
dict.ColorGMSHealthy = 'colorGMSHealthy' ;
dict.InSpeed = 'inSpeed' ;
dict.OutLoad = 'outLoad' ;
end % of function 'buildDictionary'