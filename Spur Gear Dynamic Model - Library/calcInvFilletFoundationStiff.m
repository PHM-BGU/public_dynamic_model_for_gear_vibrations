function invKf = calcInvFilletFoundationStiff(alpha, Wx, E, X, Y, dedendumR, holeR, coeffMtx, defectInfo)
%{
% Description:
% This function calculates the inverse stiffness of the residual stress
% induced by the fillet foundation. The inverse stiffness is calculated
% empirically by diving the displacement of the tooth by the contact force.
% For further information about the calculations, please see manual.
% =====
% Inputs:
% * alpha - [rad], pressure angle, may be a vector.
% * Wx - tooth width [mm] - vector along tooth axis (x axis)
% * E - [MPa], Young modulus of elasticity
% * [X,Y] - [mm], Cartesian coordinates of the tooth profile.
% * dedendumR - [mm], the radius of the dedendum circle.
% * holeR - [mm], the radius of the hole.
% * coeffMtx - coefficient matrix for the empirical calculations.
% * defectInfo - a structure with relevant information about the fault.
% =====
% Outputs:
% invKf - [N/mm]^-1, The inverse of the fillet foundation stiffness.
% =====
% In-Function Variables:
% * F - Equivalent (unit) contact force.
% * h - The ratio between the dedendum radius and the hole radius.
% * u - [mm], the radial distance between the dedendum circle and the intersection
%   point of the pressure line with the center line of the tooth profile.
% * theta_f - [rad], The angle between the center line of the tooth profile
%   and the intersection point of the profile on the dedendum circle.
%   Recall that the first index of [X,Y] is at the root.
% * Sf - [mm], Length of the arc between the edges of the tooth profile
%   on the dedendum circle.
% * polyVars - the independent varaibles in the polynoms.
% * valsLMPQ - the values for [L,M,P,Q], respectively.
% * delta_f - The displacement induced by fillet foundation.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Handle recursively a special case of ToothFaceFault through %%
% The through tooth face fault affects the stiffness of the healthy tooth
% of the other wheel due to the varying pressure angle and contact points.
if strcmp(defectInfo.status,'Healthy') && ...
        strcmp(defectInfo.statusOtherWheel, 'ThroughFaceFault')
    defectInfo.statusOtherWheel = 'Healthy' ;
    defectInfo.status = 'ThroughFaceFault' ;
    defectInfo.Y = Y ;
    invKf = calcInvFilletFoundationStiff(alpha, Wx, E, X, Y, dedendumR, holeR, coeffMtx, defectInfo) ;
    return
end % of if

%% Modify [X,Y] in case of a through ToothFaceFault %%
try
    [startInd, midInd, endInd] = ...
        deal(defectInfo.faultInds.start,...
        defectInfo.faultInds.mid, defectInfo.faultInds.end) ;
    if endInd > startInd
        start2midInds = startInd:midInd ;
        mid2endInds = (midInd+1):endInd ;
    else
        start2midInds = midInd:startInd ;
        mid2endInds = endInd:(midInd-1) ;
    end
    [X(start2midInds), Y(start2midInds)] = deal(X(startInd), Y(startInd)) ;
    [X(mid2endInds), Y(mid2endInds)] = deal(X(endInd), Y(endInd)) ;
end % of try

%% Define geometrical parameters for the empirical calculations %%
F = 1 ; % apply a unit force
h = dedendumR / holeR ;
u = X - dedendumR - Y.*tan(alpha) ;
theta_f = atan(Y(1)/X(1)) ;
Sf = 2 * theta_f * dedendumR ; % multiplying by 2 due to symmetry.

%% Calculating delta_f empirically %%
polyVars = [1/theta_f^2, h^2, h/theta_f, 1/theta_f, h, 1].' ; % the independent varaibles
valsLMPQ = coeffMtx * polyVars ;
[L, M, P, Q] = deal(valsLMPQ(1), valsLMPQ(2), valsLMPQ(3), valsLMPQ(4)) ;
delta_f = (F.*cos(alpha).^2./(E*Wx)).*(...
    (L.*(u./Sf).^2) + M.*(u./Sf) + P.*(1+Q.*tan(alpha).^2)) ;

%% Express the inverse stiffness with the empirical term %%
invKf = delta_f/F ;

end % of function "calcInvFilletFoundationStiff"