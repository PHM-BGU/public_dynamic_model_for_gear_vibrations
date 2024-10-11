function [wheelDefectInfo, otherWheelDefectInfo] = ...
    genDefectInfoStruct(wheelInfo, wheelDefectInfo, otherWheelInfo, otherWheelDefectInfo)
%{
% Description:
% This function generates the parameters for each fault, according
% to the fault's requirements.
% =====
% Inputs:
% * other/wheelInfo - a structure with all the information of the healthy tooth.
% * other/wheelDefectInfo - a structure with the fault type and dimesions.
% =====
% Outputs:
% * other/wheelDefectInfo - updated with the relevant information for the simulation.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Input check %%
if strcmp(wheelInfo.wheel, 'in')
    error('The model does not deal with faults on the driving wheel')
end % of if

%% Store vital information in wheelDefectInfo %%
wheelDefectInfo.otherWheel = otherWheelInfo ;
wheelDefectInfo.initContInd = wheelInfo.initContInd ;

%% Handle each defect case %%
switch wheelDefectInfo.status
    case 'Healthy'
        return
    case 'ToothBreakage'
        wheelDefectInfo = genToothBreakageProfile(wheelInfo, wheelDefectInfo) ;
    case {'PartialFaceFault', 'ThroughFaceFault'}
        wheelDefectInfo = genToothFaceFaultProfile(wheelInfo, wheelDefectInfo) ;
    case 'ToothDestruction'
        wheelDefectInfo = genToothDestructionProfile(wheelInfo, wheelDefectInfo) ;
end % of switch-case

%% Generate the defect info structure of the other undamaged wheel %%
otherWheelDefectInfo = wheelDefectInfo.otherWheel ;
otherWheelDefectInfo.status = 'Healthy' ;
otherWheelDefectInfo.statusOtherWheel = wheelDefectInfo.status ;
otherWheelDefectInfo.zOtherWheel = wheelInfo.z ;
otherWheelDefectInfo.otherWheel = wheelDefectInfo ;
otherWheelDefectInfo.otherWheel.otherWheel = [] ; % remove unnecessary fields

end % of function 'genDefectInfoStruct'