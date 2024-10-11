function defectInfo = genToothBreakageProfile(wheelInfo, defectInfo)
%{
% Description:
% This function calculates the profile of a broken tooth face,
% given the healthy tooth profile and the defect parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      x ^                %
%                      __|________        %
%          Fault->   /   |        |       %
%                  /     |        |       %
%                /       |        |       %
%               |        |        |       %
%               |        |        |       %
%%     ------------------|----------- >z %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====
% Inputs:
% * wheelInfo.X - [mm], X coordinates of the tooth.
% * wheelInfo.W - [mm], tooth width.
% * wheelInfo.initContInd - initial contact index.
% * wheelInfo.epsilon - contact ratio.
% * defectInfo - defect parameters, relevant fields for this function:
%   # defectInfo.zOtherWheel - the number of teeth on the other wheel.
%   # defectInfo.dimensions.x/z - [mm], the height (x) and width (z) of
%   the plain (line in paractice) that chipps the tooth.
%   # defectInfo.dimensions.tipLoss - [mm], tooth height removed from the
%   tip. A pivotal parameter that does not participate in this particulat
%   function but will be encountered later.
% =====
% Output:
% * defectInfo.defEntranceCycIn - an index of the initial phase (theta)
%   of the fault. Since the fault starts at the tooth tip, the initial 
%   phase is set to zero (0).
% * defectInfo.defExitCycIn - an index of the phase (theta) at the
%   final point of the fault.
% * defectInfo.Wx - [mm], a vector with the tooth width (W) along X.
% * defectInfo.defWidth - [mm], a vector with the missing width induced by the defect.
% * defectInfo.Zc - [mm], a vector with the z-coordinates of the center of area.
% * defectInfo.otherWheel.(Wx/Zc/defWidth) - [mm], same parameters
%   calculated for the other (healthy) wheel
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Load fault dimensions and basic parameters %%
X = wheelInfo.X ;
zOtherWheel = defectInfo.zOtherWheel ;
xDef = defectInfo.dimensions.x ;
zDef = defectInfo.dimensions.z ;

%% Calculate fault's entrance & exit points according to the input cycle %%
toothTipInd = length(X) ;
chipEdgeInd = find(X > X(end) - xDef, 1) - 1 ; % index to the chipping edge
initContInd = wheelInfo.initContInd ;

defectInfo.defEntranceCycIn = 0 ; % the breakage starts at the tooth tip
defectInfo.defExitCycIn = ...
    (toothTipInd - chipEdgeInd) / ...
    (toothTipInd - initContInd) * ...
    (wheelInfo.epsilon / zOtherWheel) ;

%% Calculate tooth width along X %%
[W0, defectInfo.otherWheel.Wx] = deal(wheelInfo.W) ;
[defWidth, defectInfo.otherWheel.defWidth] = deal(zeros(size(X))) ;
defWidth(chipEdgeInd:end) = (zDef/xDef) * (X(chipEdgeInd:end) - X(chipEdgeInd)) ;
defectInfo.defWidth = defWidth ;
defectInfo.Wx = W0 - defWidth ;
defectInfo.Zc = defWidth / 2 ;

defWidth = flip(defWidth) ;
defWidth = defWidth(1:end-initContInd+1) ;
% defWidth = resampling(defWidth, length(X)-defectInfo.otherWheel.initContInd+1) ;
defectInfo.otherWheel.defWidth(defectInfo.otherWheel.initContInd:end) = defWidth ;
defectInfo.otherWheel.Zc = 0 ;
end % of function 'genToothBreakageProfile'