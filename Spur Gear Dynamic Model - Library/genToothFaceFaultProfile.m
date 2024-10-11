function defectInfo = genToothFaceFaultProfile(wheelInfo, defectInfo)
%{
% Description:
% This function calculates the profile of a tooth with a face fault,
% given the healthy tooth profile and the defect parameters.
% Additional relevant parameters are calculated depending on the fault
% type, that is, PartialFaceFault or ThroughFaceFault.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     x ^                %
%                     __|__              %
%                    /  |  \             %
%                    \  |   \            %
%          Fault->    | |    \           %
%                  _ /  |     \          %
%              __/      |      \__       %
%%     y<---------------|-------------  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====
% Inputs:
% * [wheelInfo.X, wheelInfo.Y] - [mm], coordinates of a healthy profile.
% * wheelInfo.alpha - [rad], nominal pressure angle.
% * wheelInfo.epsilon - contact ratio.
% * wheelInfo.initContInd - index to the initial contact point.
% * wheelInfo.W - [mm], nominal tooth width.
% * defectInfo - defect parameters, relevant fields for this function:
%   # defectInfo.dimensions - the parameters of the defecting circle, see 
%   'calcDefectedFaceProfile' or assist user manual for more information. 
%   # defectInfo.zInWheel - number of teeth on the driving wheel, necessary
%   for calculations in the cycls domain.
% =====
% Output:
% * defectInfo.Y - [mm], the defected tooth profile.
% * defectInfo.faultInds.start/mid/end - indices in X to the fault edges.
% * defectInfo.defLen - [mm], the length of the line connecting the fault
% edges, can be vital for data labeling.
% * defectInfo.faultCycs - cycle points of the input shaft corresponding to
% the fault edges.
% * defectInfo.Wx - [mm], variation of the tooth width along X.
% * defectInfo.defWidth - [mm], the missing width induced by the defect.
% * defectInfo.geometry - see calcToothFaceFaultGeometry or assist user
% manual for more information.
% * defectInfo.otherWheel - all the fields above but calculated for the
% other (healthy) wheel.
% * defectInfo.h - [mm], deviation from involute, relevant for calculating
% the mesh forces induced by fault displacement.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Settings of the healthy tooth profile %%
[X, Y, initContInd] = deal(wheelInfo.X, wheelInfo.Y, wheelInfo.initContInd) ;

%% Calculate defY and fault indices in X %%
[defY, defStartInd, defEndInd] = ...
    calcDefectedFaceProfile(X, Y, initContInd, defectInfo.dimensions) ;

[defectInfo.Y, defectInfo.faultInds.start, defectInfo.faultInds.end] = ...
    deal(defY, defStartInd, defEndInd) ;
defectInfo.faultInds.mid = floor(0.5*(defStartInd+defEndInd)) ;
defectInfo.defLen = ...
    norm([X(defStartInd) Y(defStartInd)] - [X(defEndInd) Y(defEndInd)]) ;
defectInfo.defCoeff = defectInfo.defLen / norm([X(initContInd) Y(initContInd)] - [X(end) Y(end)]) ;

%% Calculate the input cycle points at the start/end of the fault %%
defectInfo.faultCycs.start = (wheelInfo.epsilon/defectInfo.zInWheel) * ...
    (X(defEndInd) - X(end))/(X(initContInd)-X(end)) ; % [cycle]
defectInfo.faultCycs.end = (wheelInfo.epsilon/defectInfo.zInWheel)* ...
    (X(defStartInd) - X(end))/(X(initContInd)-X(end)) ; % [cycle]

%% Calculate fault indices to the other (healthy) wheel %%
lenX = length(X) ;
defectInfo.initContInd = initContInd ;
indLabels = {'start', 'mid', 'end'} ;
for ii = 1:length(indLabels)
    defectInfo.otherWheel.faultInds.(indLabels{ii}) = ...
        lenX + initContInd - defectInfo.faultInds.(indLabels{ii}) ;
end % of for ii

%% Calculate relevant parameters for partial/through face fault %%
switch defectInfo.status
    case 'PartialFaceFault'
        %% (PartialFaceFault) Calculate variation of tooth width and defect width along X %%
        [W, defectInfo.otherWheel.Wx] = deal(wheelInfo.W) ;
        defW = defectInfo.dimensions.W ;
        [defWidth, defectInfo.otherWheel.defWidth] = deal(zeros(size(X))) ;
        defWidth(defStartInd:defEndInd) = defW ;
        defectInfo.defWidth = defWidth ;
        defectInfo.Wx = W - defWidth ;
        
        %% (PartialFaceFault) Repeat on the former block for the other (healthy) wheel %%
        defWidth = flip(defWidth) ;
        defWidth = defWidth(1:end-initContInd+1) ;
        defectInfo.otherWheel.defWidth(defectInfo.otherWheel.initContInd:end) = defWidth ;
        
        %% (PartialFaceFault) Calculate cross-section properties %%
        [defectInfo.geometry.A, ...
            defectInfo.geometry.Yc, defectInfo.geometry.Zc, ...
            defectInfo.geometry.Iy, defectInfo.geometry.Iz, ...
            defectInfo.geometry.Iyz] = ...
            calcToothFaceFaultGeometry(Y, defY, W, defW) ;
        
    case 'ThroughFaceFault'
        %% (ThroughFaceFault) Calculate deviation from involute %%
        hWheel = calcDeviationFromCurve(X, Y, defStartInd, defEndInd) ;
        otherWheel = defectInfo.otherWheel ;
        hOtherWheel = calcDeviationFromCurve(otherWheel.X, otherWheel.Y, ...
            otherWheel.faultInds.start, otherWheel.faultInds.end) ;
        defectInfo.h = hOtherWheel + hWheel ;
        
        %% (ThroughFaceFault) Calculate variation of alpha along X %%
        defectInfo = calcAlphaX(wheelInfo.alpha, lenX, defectInfo) ;
end % of switch defectInfo.status
end % of function 'genToothFaceFaultProfile'
%%
function [defY, defStartInd, defEndInd] = calcDefectedFaceProfile(X, Y, initContInd, circParam)
%{
% Description:
% This function calculates the profile of a tooth with a face fault,
% given the healthy tooth profile and the defect parameters.
% =====
% Inputs:
% * [X,Y] - [mm], Cartesian coordinates of the healthy tooth face.
% * initContInd - initial contact index of the involute profile.
% * circParam - a structure with parameters of the defect circle.
% =====
% Outputs:
% * defY - [mm], defected tooth face profile.
% * [defStartInd, defEndInd] - fault edge indices.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
%% Calculate the radius and circle center coordinates %%
switch circParam.method
    case 'byPercentage'
        lenCont = length(X)-initContInd ;
        defStartInd = initContInd + round(circParam.defStartCoeff * lenCont) ;
        defEndInd = initContInd + round(circParam.defEndCoeff * lenCont) ;
        [xCen, yCen, R] = findCircCenAndR( ...
            [X(defStartInd) ; Y(defStartInd)], [X(defEndInd) ; Y(defEndInd)], ...
            circParam.depthCoeff) ;
    case 'byCircCoordinates'
        [xCen, yCen, R] = deal(circParam.circCen.x, circParam.circCen.y, circParam.circR) ;
    otherwise
        error('ToothFaceFault is under-defined')
end % of switch-case

%% Calculate the Cartesian coordinates of the lower semicircle %%
theta = linspace(pi, 2*pi, circParam.circNPts)' ;
xCirc = xCen + R*cos(theta) ;
yCirc = yCen + R*sin(theta) ;

%% Calculate fault edge indices %%
[xCircLowerBound, xCircUpperBound] = deal(xCen-R, xCen+R) ;
[~, circStartInd] = min(abs(X - xCircLowerBound)) ;
circStartInd = max(circStartInd, initContInd) ;
[~, circEndInd] = min(abs(X - xCircUpperBound)) ;
[xInvlt, yInvlt] = deal(X(circStartInd:circEndInd), Y(circStartInd:circEndInd)) ;
yInvlt = interp1(xInvlt, yInvlt, xCirc, 'linear', 'extrap') ;
[~, distInds] = sort(abs(yInvlt - yCirc)) ; % the first two inds are the fault edges.
distInds(diff(distInds(1:2)) == 1) = [] ; % avoid consecutive points
faultEdges = yInvlt(sort(distInds(1:2))) ;

%% Calculate the damaged tooth profile (defY) %%
defStartInd = find(Y>faultEdges(1), 1, 'last') + 1 ; % (+1) to assure index within the circle
defEndInd = find(Y>faultEdges(2), 1, 'last')  - 1 ; % (-1) to assure index within the circle
defY = Y ;
defY(defStartInd:defEndInd) = interp1(xCirc, yCirc, X(defStartInd:defEndInd), 'spline') ;
end % of function 'calcDefectedFaceProfile'
%%
function [xCen, yCen, radius] = findCircCenAndR(pt1, pt2, depthCoeff)
%{
% Description:
% This function calculates the radius and the center of a circle based on
% two points on its circumference a depth coefficient (between 0 to 1).
% =====
% Inputs:
% * pt1/pt2 - point on the circle's circumference [x;y].
% * depthCoeff - [0,1] - please assist user manual for details.
% =====
% Output:
% * [xCen, yCen] - [x,y] cooridnates of the circle's center.
% * radius - circle's radius.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
%% Calculate chord length and circle radius %%
chord = norm(pt1 - pt2) ;
radius = chord/4 * (depthCoeff + 1/depthCoeff) ;

%% Calculate the chord's perpendicular bisector (magnitude and direction) %%
bisectorMag = sqrt(radius^2 - (chord/2)^2) ;
bisectorDir = abs(flip(pt2 - pt1) .* [-1 ; 1] / chord) ; % abs ensures 1st quadrant

%% Calculate the circle center (upper point among two options) %%
midchord = (pt1 + pt2)/2 ;
upperCen = midchord + bisectorMag*bisectorDir ;
[xCen, yCen] = deal(upperCen(1), upperCen(2)) ;
end % of function 'findCircCenAndR'
%%
function [A, Yc, Zc, Iy, Iz, Iyz] = calcToothFaceFaultGeometry(Y, defY, W, defW, z0)
%{
% Description:
% This function calculates the cross-sectional properties of a tooth
% with a partial tooth face fault. In the zy plane, the cross-section has
% the shape of a rectangle of a healthy tooth, with a subtraction of a
% smaller rectangle from the top edge to represent the fault.
%
% %%%% Illustration of the cross-section in the y-z plane: %%%%
%
%                       y ^
%      z0                 |
%     ****         *******|*********************
%     *  *  defW   *      |          ^         *
%     *  ***********      |          |         *
%     *       ^           |          | Y       *
%     *       | defY      |          |         *
%   --------------------------------------------------->z
%     *                   |          |         *
%     *                   |          | Y       *
%     *                   |          |         *
%     *                   |          v         *
%     ********************|*********************
%
%     <------------------ W ------------------->
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====
% Inputs: (see illustration)
% * Y - [mm], Height of the healthy tooth face.
% * defY - [mm], Height of the faulted tooth face.
% * W - [mm], Width of the healthy tooth.
% * defW - [mm], Defect width.
% * z0 - [mm], Reference z-coordinate of the fault edge (positive value).
%   If this input is missing, the function will set it defaultly to zero.
%   The user is also able to set it as 'middle' and the offset will be
%   calculated within the function.
% * z0 - [mm], Reference z-coordinate of the fault edge (positive value).
%   If this input is missing, the function will default it to zero.
%   Alternatively, you can set it as 'middle,' and the function will
%   calculate the offset automatically to position the fault edge at the
%   midpoint of the tooth width.
% =====
% Outputs:
% * A - [mm^2], Cross-sectional area.
% * Yc - [mm], Y-coordinate of the center of area.
% * Zc - [mm], Z-coordinate of the center of area.
% * Iy, Iz, Iyz - [mm^4], Area moments of inertia with respect to a
%   parallel coordinate system in the center of area.
% =====
% In-Function Variables:
% * h$ - Property of the healthy tooth area.
% * def$ - Property of the defected area.
% =====
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Input check %%
if nargin == 4
    z0 = 0 ;
elseif strcmp(z0, 'middle')
    z0 = (W-defW)/2 ;
end % of if

%% Calculate cross-section area %%
hA = W.*(2*Y) ;
defA = defW.*(Y - defY) ;
A = hA - defA ;

%% Calculate Y-coordinate of the center of area %%
hYc = 0 ;
defYc = 0.5*(Y + defY) ;
Yc = (hYc.*hA - defYc.*defA) ./ A ;

%% Calculate Z-coordinate of the center of area %%
hZc = 0 ;
defZc = z0 - 0.5*(W - defW) ;
Zc = (hZc.*hA - defZc.*defA) ./ A ;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Clarification: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The area moments of inertia are calculated in the former steps %%%%%
%%%%% using the principle of superposition: I = hI - defI , %%%%%%%%%%%%%%
%%%%%     and Steiner's theorem: I = Ic + delta1*delta2*A   . %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate area moment of inertia about the Y-axis %%
hIy = (1/12*W.^3.*(2*Y)) + (hZc - Zc).^2.*hA ;
defIy = (1/12*defW.^3.*(Y - defY)) + (defZc - Zc).^2.*defA ;
Iy = hIy - defIy ;

%% Calculate area moment of inertia about the Z-axis %%
hIz = (1/12*W.*(2*Y).^3) + (hYc - Yc).^2.*hA ;
defIz = (1/12*defW.*(Y - defY).^3) + (defYc - Yc).^2.*defA ;
Iz = hIz - defIz ;

%% Calculate product moment of area %%
hIyz = 0 + (hZc - Zc).*(hYc - Yc).*hA ;
defIyz = 0 + (defZc - Zc).*(defYc - Yc).*defA ;
Iyz = hIyz - defIyz ;

end % of function 'calcToothFaceFaultGeometry'
%%
function h = calcDeviationFromCurve(X, Y, ind1, ind2)
%{
% Description:
% This function calculates the maximal deviation h of a curve from a line
% stretched between its boundaries indexes [ind1,ind2]. The coordinate
% system is transformed to a system wherein the stretched line is
% horizontal in aspect to the original y axis. It is assumed that after
% rotation, the minimum y value equals to the y value of the stretched line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     [X1,Y1] o**************************o [X2,Y2]     %%
%%             *           ^             *              %%
%%              *          | h          *               %%
%%               *         |          *                 %%
%%                 *       v       *                    %%
%%                     *********                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====
% Inputs:
% * [X,Y] - [mm], original coordinates of the curve.
% * [ind1, ind2] - indexes of the curve's boundaries.
% =====
% Output:
% * h - [mm], maximal deviation from curve.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
if ind1>ind2
    [ind1,ind2] = deal(ind2,ind1) ;
end % of if
phi = atan2(abs(Y(ind1)-Y(ind2)), abs(X(ind1)-X(ind2))) ; % rotation angle
rotMtx = [cos(phi) -sin(phi) ; sin(phi) cos(phi)] ; % rotation matrix
yRot = [0,1]*rotMtx*[X(ind1:ind2), Y(ind1:ind2)]' ; % transformed y points.
h = max(yRot) - min(yRot) ;
end % of function 'calcDeviationFromCurve'.
%%
function defectInfo = calcAlphaX(alphaNom, lenX, defectInfo)
%{
% Description:
% This function computes the variable pressure angle in case of a throguh
% tooth face fault. It is designed to model the scenario where the healthy
% wheel's tooth deviates from the nominal pressure line and instead rotates
% around the fault edges, with a transition occurring in the middle between these edges.
% =====
% Inputs:
% * alphaNom - [rad], nominal pressure angle.
% * lenX - length of the tooth profile (same length of the output).
% * defectInfo - relevant fields for this function (surprisingly, most of
% the required info belongs to the healthy tooth rather than the damaged).
%   # defectInfo.faultInds.start/end - start and end indices of the fault.
%   # defectInfo.otherWheel - information of the other (healthy) wheel:
%       o otherWheel.R.base - [mm], base radius.
%       o otherWheel.faultInds.start/end - start and end indices of the
%         fault on the other (healthy) wheel.
%       o otherWheel.X/Y - [mm], coordinates of the other (healthy) wheel.
% =====
% Output:
% * defectInfo.alphaX  - [rad], alpha(X_wheel), with the same length as X.
% * defectInfo.otherWheel.alphaX - [rad], pressure angle of the other wheel.
% * defectInfo.defR - [mm], a structure with the relevant radii:
%   # startDefR - The radius of the other wheel at the fault's start.
%   # midDefR - The radius of the other wheel at the middle of the fault.
%   # endDefR - The radius of the other wheel at the fault's exit.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Calculate radii of the fault on the other (healthy) wheel %%
otherWheel = defectInfo.otherWheel ;
faultIndsOther = defectInfo.otherWheel.faultInds ;
Rb = otherWheel.R.base ;
startDefR = norm([otherWheel.X(faultIndsOther.start), otherWheel.Y(faultIndsOther.start)]) ;
endDefR = norm([otherWheel.X(faultIndsOther.end), otherWheel.Y(faultIndsOther.end)]) ;
midDefR = norm([otherWheel.X(faultIndsOther.mid), otherWheel.Y(faultIndsOther.mid)]) ;
[defectInfo.defR.start, defectInfo.defR.mid, defectInfo.defR.end] = deal(startDefR, midDefR, endDefR) ;

%% Build vectors for two ranges of indices (transition at the midpoint) %%
faultInds = (defectInfo.faultInds.start:defectInfo.faultInds.end)' ;
lenRx = ceil(length(faultInds)/2) ;
mid2startRx = linspace(startDefR, midDefR, lenRx) ;
end2midRx = linspace(midDefR, endDefR, lenRx) ;

%% Calculate the variable pressure angle %%
alphaX = alphaNom*ones(lenX, 1) ;
for ii = 1:lenRx
    alphaX(faultInds(ii)) = alphaNom + ...
        (sqrt(startDefR^2-Rb^2)-sqrt(mid2startRx(ii)^2-Rb^2))/Rb ;
    alphaX(faultInds(ii+lenRx-1))= alphaNom + ...
        (sqrt(endDefR^2-Rb^2)-sqrt(end2midRx(ii)^2-Rb^2))/Rb ;
end % of for

%% Update both the angle of the wheel and other wheel %%
defectInfo.alphaX = alphaX ;
alphaXOther = alphaNom * ones(lenX, 1) ;
faultIndsSrt = sort([faultIndsOther.start, faultIndsOther.end]) ;
alphaXOther(faultIndsSrt(1):faultIndsSrt(2)) = resampling( ...
    flip(alphaX(faultInds)), length(faultIndsSrt(1):faultIndsSrt(2)), 'PCHIP') ;
defectInfo.otherWheel.alphaX = alphaXOther ;
end % of function 'calcAlphaX'
%%