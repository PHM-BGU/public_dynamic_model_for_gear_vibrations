function [inWheelInfo, outWheelInfo] = calcContLineProp(inWheelInfo, outWheelInfo)
%{
% Description:
% This function calculates the initial contact index for each wheel and the
% contact ratio for spur gears. In addition, this function matches the
% dimensions of the tooth profile of the input and output wheels in
% the contact length, between the initial contact point to the tip.
% =====
% Inputs:
% * (in/out)WheelInfo - the relevant information about the wheels:
%   # (in/out)WheelInfo.z - number of teeth.
%   # (in/out)WheelInfo.alpha - [rad], pressure angle.
%   # (in/out)WheelInfo.R - [mm], a struct with all the radii.
%   # (in/out)WheelInfo.(X/Y) - [mm], tooth profile coordinates.
% =====
% Outputs:
% * (in/out)WheelInfo - updated with three pivotal fields:
%   # initContInd - initial contact index along the involute profile.
%   # R.initCont - the radius at the initial contact point.
%   # epsilon - contact ratio.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
%% Find the initial contact point for each wheel %%
[inInitContInd, inWheelInfo.R.initCont] = findInitContInd(inWheelInfo, outWheelInfo) ;
[outInitContInd, outWheelInfo.R.initCont] = findInitContInd(outWheelInfo, inWheelInfo) ;

%% Match the dimensions of the contact points (from initial contact to tip) %%
inXY = [inWheelInfo.X, inWheelInfo.Y] ;
[inRootInds, inContInds] = deal(1:inInitContInd-1, inInitContInd:length(inWheelInfo.X)) ;
[lenRootInds, lenContInds] = deal(outInitContInd-1, length(inWheelInfo.X)-outInitContInd+1) ;
inRootXY = resampling(inXY(inRootInds,:), lenRootInds, 'PCHIP') ;
inContXY = resampling(inXY(inContInds,:), lenContInds, 'PCHIP') ;
inWheelInfo.X = [inRootXY(:,1) ; inContXY(:,1)] ;
inWheelInfo.Y = [inRootXY(:,2) ; inContXY(:,2)] ;
[inWheelInfo.initContInd, outWheelInfo.initContInd] = deal(outInitContInd) ;

%% Calculate the contact ratio ε %%
epsilon = calcContactRatio(inWheelInfo.z, outWheelInfo.z, inWheelInfo.alpha) ;
[inWheelInfo.epsilon, outWheelInfo.epsilon] = deal(epsilon) ;
end % of function 'calcContLineProp'
%%
function [initContInd, r] = findInitContInd(wheel1, wheel2)
%{
% Description:
% This function calculates the initial contact index based on the
% geometric initial contact point on the involute near the base circle.
% =====
% Inputs:
% * wheel1/wheel2 - a structure with all the geometric parameters:
%   # wheel1/wheel2.alpha - [rad], pressure angle.
%   # wheel1/wheel2.R - [mm], a struct with all the radii.
%   # wheel1/wheel2.(X/Y) - [mm], tooth profile coordinates.
% =====
% Outputs:
% * initContInd - an index of the initial contact point of wheel1.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
[alpha, Ra2, Rb2, Rp2, Rb1, Rp1] = deal(wheel1.alpha, wheel2.R.addendum, ...
    wheel2.R.base, wheel2.R.pitch, wheel1.R.base, wheel1.R.pitch) ;
r = sqrt(...
    (sqrt(Ra2^2-Rb2^2)-(Rp1+Rp2)*sin(alpha))^2 + ...
    Rb1^2 ) ;
rInvlt = sqrt(wheel1.X.^2 + wheel1.Y.^2) ;
[~,initContInd] = min(abs(rInvlt-r)) ; % argmin
end % of function 'findInitContInd'
%%
function epsilon = calcContactRatio(z1, z2, alpha)
%{
% Description:
% This function calculates the contact ratio of two gears.
% =====
% Inputs:
% z1 - number of teeth on wheel 1 (does not matter if driving/driven)
% z2 - number of teeth on wheel 2
% alpha - [rad], pressure angle.
% =====
% Outputs:
% epsilon - contact ratio
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
element1 = sqrt( (z1+2)^2 - (z1*cos(alpha))^2 ) ;
element2 = sqrt( (z2+2)^2 - (z2*cos(alpha))^2 ) ;
element3 = (z1+z2) * sin(alpha) ;
epsilon = (element1 + element2 - element3) / (2*pi * cos(alpha)) ;
end % of function 'calcContactRatio'
%%