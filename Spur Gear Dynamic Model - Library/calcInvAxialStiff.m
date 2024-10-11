function invKa = calcInvAxialStiff(alpha, W, E, X, Y, defectInfo)
%{
% Description:
% This function calculates the inverse stiffness derived from the
% potential strain energy of the tooth due to compressive axial stres.
% The energy is calculated by integration using the cumtrapz function.
% =====
% Inputs:
% * alpha - [rad], pressure angle, may be a vector alpha(X).
% * W - [mm], tooth width.
% * E - [MPa], Young modulus of elasticity.
% * [X,Y] - [mm], Cartesian coordinates of a healthy tooth profile.
% * defectInfo - a structure with relevant information about the fault.
% =====
% Outputs:
% * invKa - [N/mm]^-1, A vector of the inverse axial stiffness.
% =====
% Significant In-Function Variables:
% * defY - [mm], The defected tooth profile
% * A - [mm^2], cross-section area (a rectangle in a healthy status).
% * Ua - potential energy due to axial stresses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Handle recursively a special case of TroughFaceFault %%
%%% This step considers the case of a healthy tooth, while the
%%% tooth of the other wheel is damaged.
if strcmp(defectInfo.status,'Healthy') && strcmp(defectInfo.statusOtherWheel,'ThroughFaceFault')
    defectInfo.statusOtherWheel = 'Healthy' ;
    defectInfo.status = 'ThroughFaceFault' ;
    defectInfo.Y = Y ;
    invKa = calcInvAxialStiff(alpha, W, E, X, Y, defectInfo) ;
    return
end

%% Apply a unit force %%
F = 1 ; % unit contact force
Fa = F * sin(alpha) ; % axial component of the contact force

%% Calculate the potential energy due to axial (tension/compression) stress Ua %%
switch defectInfo.status
    case 'Healthy'
        A = W.*(2*Y) ;
        Ua = calcUa(X, E, A, Fa) ;
        
    case 'ToothBreakage'
        A = defectInfo.Wx.*(2*Y) ;
        Ua = calcUa(X, E, A, Fa) ;
        
    case 'PartialFaceFault'
        A = defectInfo.geometry.A ;
        Ua = calcUa(X, E, A, Fa) ;
        
    case 'ThroughFaceFault'
        defY = defectInfo.Y ;
        A = W.*(Y+defY) ;
        Ua = calcUa(X, E, A, Fa) ;
        Ua = correctUForThroughFaceFault(Ua, defectInfo.faultInds, Fa) ;
        
    case 'ToothDestruction'
        defY = defectInfo.Y ;
        A = W.*(Y+defY) ;
        Ua = calcUa(X, E, A, Fa) ;
        
end % of switch defectInfo.status

%% Express the inverse stiffness with the potential energy %%
invKa = (2*Ua)/(F^2) ;

end % of function "calcInvAxialStiff"
%%
function Ua = calcUa(X, E, A, Fa)
%{
% Description:
% This function computes the potential strain energy resulting from
% compressive axial stress according to beam theory.
% Improtant note: The axial force Fa may vary along X but remains
% X-independent within the cumulative integral. Consequently, the term
% involving Fa is extracted from the integral and replaced by scalar multiplication.
% =====
% Inputs & Outputs:
% * please see the description in the parent function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
Ua = (Fa.^2).*cumtrapz(X, 1./(2*E*A)) ;
end % of function 'calcUa'