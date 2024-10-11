function [xTop, yTop, R, relEdgeInd, xBot, yBot, xAdd, yAdd] = ...
    genHealthyToothProfile(module, z, alpha_0, toothModification, NPtList, thetaInit)
%{
% Description:
% This function takes gear wheel parameters (metric module, number of teeth,
% and pressure angle) and calculates the tooth profile in Cartesian
% coordinates along with the characteristic radii. Optional tooth profile
% modifications, such as tip relief and crowning, can be incorporated upon request.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          x ^                    %
% addendum (tip) ->        __|__                  %
%                         /  |  \                 %
% operational            /   |   \      inactive  %
% face (Top) ->         /    |    \      <- face  %
%                      /     |     \       (Bot)  %
% dedendum (root)-> __/      |      \__           %
%%          y<---------------|-------------      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====
% Inputs:
% * module - [mm], metric module of the gear.
% * z - number of teeth.
% * alpha_0 - [rad], pressure angle.
% * NPtList - a list with the resolution of the profile
% * thetaInit - the initial angle of the tooth (good for plottings).
% =====
% Outputs:
% * [xTop, yTop] - The Cartesian horizontal and vertical coordinates
%   of the profile of the operational tooth face (see illustration above),
%   starting at the dedendum circle, ending at the addendum circle. [mm]
% * [xBot, yBot] - The Cartesian horizontal and vertical coordinates
%   of the profile of the inactive tooth face (see illustration above),
%   starting at the dedendum circle, ending at the addendum circle. [mm]
% * R - A list structure holding the characteristic radii of the gear:
%     # R.pitch - pitch circle radius [mm]
%     # R.addendum - addendum circle radius [mm]
%     # R.dedendum - dedendum circle radius [mm]
%     # R.base - base circle radius [mm]
% * [xAdd, yAdd] - [mm], Tooth tip coordinates, i.e., arc on the addendum circle
% * relEdgeInd - an index corresponding to the edge of the relieved tip.
% =====
% In-Function Variables:
% * pitchDiam - Pitch circle diameter [mm]
% * circPitch - Circular pitch, i.e., the arc length of one tooth on the
%   pitch circle. Half of the circular pitch is the tooth thickness,
%   while the other half consists of the rest of the arc of the tooth pattern
% * s0 - Tooth thickness = arc length on the pitch circle from face to face.
% * theta_0 - The angle of the tooth thickness on the pitch circle = pi/2z.
% * invAlpha_0 - The involute function of the pressure angle alpha_0. [rad]
% * alpha_P - A vector of the pressure angle of an arbitrary point P on the involute curve.
% * invAlphaP - The involute function of the pressure angles vector alpha_P. [rad]
% * [thetaPInvltTop/Bot, RPInvlt] - Polar coordinates of the top/bottom involute profile curves.
% * [thetaAddTop, thetaAddBot] - Upper & Lower bounds of the addendum tip curve [rad]
% * thetaAdd - A vector of the angle of the addendum tip curve [rad]
% * psiFillet - The angle of the fillet curve, from the dedendum to the base [rad]
% * rhoFillet - The radius of the fillet cut.
% * gamma - A vector runnin over psiFillet.
% * [thetaFilletTop/Bot, RFillet] - Polar coordinates of the top/bottom fillet curves.
% * [xTopTemp / yTopTemp / xBotTemp / yBotTemp] - Temporary vectors of the Cartesian cooardinates.
%   These vectors are constructed from the fillet & involute curves, and need to be smoothed and downsampled.
% * startInd - The index in which the involute curve intersects the dedendum circle. Relevant only in case Db<Dd
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
    %}
    
    %% Input check %%
    if nargin < 6
        thetaInit = 2*pi ;
        if nargin < 5
            [NPtList.invlt, NPtList.addendum, NPtList.fillet, ...
                NPtList.profile]  = deal(15e3, 5e3, 5e3, 10e3) ;
            if nargin == 3
                [toothModification.tipRelief.include, ...
                    toothModification.crowning.include] = deal(false) ;
                relEdgeInd = NPtList.profile ; % tooth tip
            end % of if
        end % of if
    end % of if
    
    %% Characteristic radii calculations %%
    pitchDiam = module * z ; % [mm], pitch circle diameter
    R.pitch = 0.50 * pitchDiam ; % [mm], pitch circle radius
    R.addendum = R.pitch + 1.00 * module ; % [mm], tip circle radius
    R.dedendum = R.pitch - 1.25 * module ; % [mm], root circle radius
    R.base = R.pitch * cos(alpha_0) ; % [mm], base circle radius
    
    %% Basic calculations along the pitch circle %%
    circPitch = pi * module ; % =pi*(pitchDiam/z) [mm]
    s0 = 0.5 * circPitch ; % [mm], reference pitch tooth_thickness
    theta_0 = s0 / pitchDiam ; % [rad], reference pitch circle = pi/(2z)
    invAlpha_0 = invltFunc(alpha_0) ;
    
    %% Top & Bottom involute curves calculations %%
    alphaP = linspace(0, acos(R.base/R.addendum), NPtList.invlt+1) ;
    invAlphaP = invltFunc(alphaP) ; % vctr
    thetaPInvltTop = thetaInit + (theta_0 + invAlpha_0 - invAlphaP) ; % vctr
    thetaPInvltBot = thetaInit - (theta_0 + invAlpha_0 - invAlphaP) ; % vctr
    RPInvlt = R.base ./ cos(alphaP) ; % vctr
    [xInvltTop, yInvltTop] = pol2cart(thetaPInvltTop, RPInvlt) ;
    [xInvltBot, yInvltBot] = pol2cart(thetaPInvltBot, RPInvlt) ;
    
    %% Addendum tip curve calculations %%
    thetaAddTop = thetaInit + (theta_0 + invAlpha_0 - invAlphaP(end)) ; % upper bound
    thetaAddBot = thetaInit - (theta_0 + invAlpha_0 - invAlphaP(end)) ; % lower bound
    thetaAdd = linspace(thetaAddTop, thetaAddBot, NPtList.addendum+1)' ; % vctr
    [xAdd, yAdd] = pol2cart(thetaAdd, R.addendum) ;
    
    %% Top & Bottom fillet curves calculations %%
    psiFillet = theta_0 - invAlpha_0 ;
    rhoFillet = (R.dedendum^2 + R.base^2 -2*R.dedendum*R.base*cos(psiFillet)) / ...
        (2*R.base*cos(psiFillet) - 2*R.dedendum) ; % fillet radius
    gamma = linspace(0, psiFillet, NPtList.fillet+1) ;
    RFillet = (R.dedendum+rhoFillet)*cos(gamma) - ...
        sqrt( ((R.dedendum+rhoFillet)*cos(gamma)).^2 - ...
        R.dedendum*(R.dedendum+2*rhoFillet) ) ;
    thetaFilletTop = thetaInit + (theta_0 + invAlpha_0 + psiFillet - gamma) ;
    thetaFilletBot = thetaInit - (theta_0 + invAlpha_0 + psiFillet - gamma) ;
    [xFilletTop, yFilletTop] = pol2cart(thetaFilletTop, RFillet) ;
    [xFilletBot, yFilletBot] = pol2cart(thetaFilletBot, RFillet) ;
    
    %% Building the tooth profile %%
    if R.dedendum < R.base
        xTopTemp = [xFilletTop(1:end-1) , xInvltTop] ;
        yTopTemp = [yFilletTop(1:end-1) , yInvltTop] ;
        xBotTemp = [xFilletBot(1:end-1) , xInvltBot] ;
        yBotTemp = [yFilletBot(1:end-1) , yInvltBot] ;
    else
        [~, startInd] = min(abs((xInvltTop.^2 + yInvltTop.^2) - R.dedendum^2)) ;
        xTopTemp = xInvltTop(startInd:end) ;
        yTopTemp = yInvltTop(startInd:end) ;
        xBotTemp = xInvltBot(startInd:end) ;
        yBotTemp = yInvltBot(startInd:end) ;
    end % of 'if-else'
    
    xTop = linspace(xTopTemp(1), xTopTemp(end), NPtList.profile).' ;
    yTop = interp1(xTopTemp, yTopTemp, xTop, 'spline') ;
    xBot = linspace(xBotTemp(1), xBotTemp(end), NPtList.profile).' ;
    yBot = interp1(xBotTemp, yBotTemp, xBot, 'spline') ;
    
    %% Add crowning (optional) %%
    if toothModification.crowning.include && ~mod(abs(thetaInit), 2*pi)
        yTop = yTop - toothModification.crowning.longAmount ;
        yBot = yBot + toothModification.crowning.longAmount ;
    end % of if
    
    %% Add tip relief (optional) %%
    if toothModification.tipRelief.include && ~mod(abs(thetaInit), 2*pi)
        [xRelief, yRelief] = deal(toothModification.tipRelief.xRelief, toothModification.tipRelief.yRelief) ;
        [yTop, relEdgeInd] = genTipRelief(xTop, yTop, xRelief, +yRelief, toothModification.tipRelief.shape) ;
        yBot = genTipRelief(xTop, yBot, xRelief, -yRelief, toothModification.tipRelief.shape) ;
        addInds = find(yAdd<=yTop(end) & yAdd>=yBot(end)) ;
        [xAdd, yAdd] = deal(xAdd(addInds), yAdd(addInds)) ;
    else
        relEdgeInd = length(xTop) ;
    end % of if
    
end % of function 'genHealthyToothProfile'
%%
function invAlpha = invltFunc(alpha)
invAlpha = tan(alpha) - alpha ; % involute angle [rad]
end % of function 'invltFunc'
%%
function [y, edgeInd] = genTipRelief(x, y, xRelief, yRelief, curveShape)
%{
% Description:
% This function generates a tip relief based on a given tooth profile, tip
% relief dimensions and curve shape.
% =====
% Inputs:
% * [x, y] - [mm], Cartesian coordinates of the tooth profile.
% * [xRelief, yRelief] - [mm], distances in the x-axis and y-axis between
% the relieved edge and the tooth tip.
% * curveShape - Options: 'linear' or 'parabolic'. A linear curve is the
% simplest, resulting in a sharp relieved edge. The parabolic curve
% guarantees smooth continuity at the relieved tip edge.
% =====
% Outputs:
% * y - The involute tooth profile updated with the tip relief curve.
% * edgeInd - an index to the tip relief edge along the tooth profile.
% =====
% In-Function Variables:
% * [xEdge, yEdge] - [mm], profile coordinates at the tip relief edge.
% * curveInds - a range of indices corresponding to the tip relief curve.
% In case of parabolic curve:
% * [xTip, yTip] - [mm], new coordinates to the relieved tooth tip.
% * slopeEdge - the slope at the tip relief edge, necessary for forcing
% continuity at the tip relief edge.
% * [A,b] - nonhomogenous linear system of 3 equations with 3 variables
% representing the coefficients of the parabolic curve.
% * coeffs - a 3x1 vector with the parabola coefficients
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
%% Find the index in x-axis to the tip relief edge %%
[~, edgeInd] = min(abs(x - x(end) + xRelief)) ;
curveInds = edgeInd:length(x) ;
[xEdge, yEdge] = deal(x(edgeInd), y(edgeInd)) ;
%% Generate the tip relief curve %%
switch curveShape
    case 'linear'
        y(curveInds) = yEdge + (- yRelief / xRelief )*(x(curveInds) - xEdge) ;
    case 'parabolic'
        [xTip, yTip] = deal(x(end), y(edgeInd)-yRelief) ;
        slopeEdge = [0 ; diff(y)./diff(x)] ;
        slopeEdge = slopeEdge(edgeInd) ;
        A = [xEdge^2, xEdge, 1 ; xTip^2, xTip, 1 ; 2*xEdge, 1, 0] ;
        b = [yEdge ; yTip ; slopeEdge] ;
        coeffs = A\b ; % solve systems of linear equations Ax=b for x
        y(curveInds) = coeffs(1)*x(curveInds).^2 + coeffs(2)*x(curveInds) + coeffs(3) ;
end % of switch-case
end % of function 'genTipRelief'
%%