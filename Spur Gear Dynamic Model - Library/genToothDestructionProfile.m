function defectInfo = genToothDestructionProfile(wheelInfo, defectInfo)
%{
% Description:
% This function calculates the profile of a destructed tooth,
% given the healthy tooth profile and the defect parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     x ^                %
%                      _|__              %
%         Fault->     / |  \             %
%                   _/  |   \            %
%                  /    |    \           %
%                 /     |     \          %
%              __/      |      \__       %
%%     y<---------------|-------------  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====
% Inputs:
% * [wheelInfo.X, wheelInfo.Y] - [mm], coordinates of a healthy profile.
% * wheelInfo.initContInd - index to the initial contact point.
% * defectInfo - defect parameters, relevant fields for this function:
%   # defectInfo.dimensions.defStartCoeff - the percentage from the contact
%     length associated with the theoretical defect starting point. It lies
%     in the open interval of (-1, 1).
%   # defectInfo.dimensions.yOffset - [mm], the y-offset from the tooth
%   profile at the defect starting point.
%   # defectInfo.dimensions.method - the desired method for the fault
%     dimensions, either 'byDefSlope' or 'byDefEndCoeff'
%   # defectInfo.dimensions.byDefEndCoeff - the percentage from the contact
%     length associated with the theoretical defect ending point. It lies
%     in the open interval of (-1, 1).
%   # defectInfo.dimensions.defSlope - [rad], the defect slope angle.

% =====
% Output:
% * defectInfo.Y - [mm], The damaged tooth profile.
% * defectInfo.faultDisplacement [mm] - The fault displacement along the x
%   axis, starting from the initial contact point to the tooth tip.
% * defectInfo.faultInds.start/end - indexes to the fault edges.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Load fault dimensions and basic parameters %%
[X, Y, initContInd] = deal(wheelInfo.X, wheelInfo.Y, wheelInfo.initContInd)  ;
defStartCoeff = defectInfo.dimensions.defStartCoeff ;
yOffset = defectInfo.dimensions.yOffset ;

if abs(defStartCoeff)>=1
    error('ToothDestruction start coefficient must be in the open interval (-1,+1)')
end % of if

%% Calculate the destructed profile %%
defStartInd = initContInd + round(defStartCoeff*(length(X)-initContInd)) ;
[xDefStart, yDefStart] = deal(X(defStartInd), Y(defStartInd)-yOffset) ;

switch defectInfo.dimensions.method
    case 'byDefSlope'
        defSlope = - abs(defectInfo.dimensions.defSlope) ; % negative
    case 'byDefEndCoeff'
        defEndCoeff = defectInfo.dimensions.defEndCoeff ;
        defEndInd = initContInd + round(defEndCoeff*(length(X)-initContInd)) ;
        [xDefEnd, yDefEnd] = deal(X(defEndInd), Y(defEndInd) - yOffset) ;
        defSlope = atan2(yDefEnd-yDefStart, xDefEnd-xDefStart) ;
end % of switch-case

defY = Y ;
defLine = tan(defSlope)*(X-xDefStart) + yDefStart ;
defStartInd = max(defStartInd, initContInd) ;
defEndInd = find(defLine<Y, 1, 'last') ;

try % if noise is desired for the defective line
    defLineNoise = uniformRand(0, defectInfo.dimensions.noiseLevel, length(defStartInd:defEndInd)) ;
catch
    defLineNoise = 0 ;
end % of try-catch

defY(defStartInd:defEndInd) = defLine(defStartInd:defEndInd) + defLineNoise ;

%% Calculate the fault displacement %%
faultDisplacement = -( Y - defY ) ;
angleTgt2Invlt = abs([0 ; atan(diff(Y)./diff(X))]) ; % [rad], denoted as φ in the user manual
faultDisplacement = faultDisplacement .* cos(defSlope) ./ cos(angleTgt2Invlt + defSlope) ;

%% Update the defect information structure %
[defectInfo.Y, defectInfo.slopeRad] = deal(defY, defSlope) ;
[defectInfo.faultInds.start, defectInfo.faultInds.end] = deal(defStartInd, defEndInd) ;
defectInfo.faultDisplacement = flip(faultDisplacement(initContInd:end)) ;

end % of function 'genToothDestructionProfile'