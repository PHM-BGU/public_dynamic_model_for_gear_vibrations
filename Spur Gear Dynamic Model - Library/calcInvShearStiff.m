function invKs = calcInvShearStiff(alpha, W, G, X, Y, defectInfo)
%{
% Description:
% This function calculates the inverse stiffness derived from the
% potential energy of the tooth due to shear stresses.
% The energy is calculated by integration using the cumtrapz function.
% The function performs non-trivial calculations that are based on
% mathematical and physical analysis, described thoroughly in the manual.
% =====
% Inputs:
% * alpha - [rad], pressure angle, may be a vector.
% * W - [mm], tooth width.
% * G - [MPa], Shear modulus.
% * A - [mm^2], cross-section area (a rectangle in a healthy status).
% * [X,Y] - [mm], Cartesian coordinates of a healthy tooth profile.
% * defectInfo - a structure with relevant information about the fault.
% =====
% Outputs:
% * invKs - [N/mm]^-1, A vector of the inverse shear stiffness.
% =====
% In-Function Variables:
% * Us - potential energy due to shear stresses.
% =====
% In-Function Variables Related to Partial ToothFaceFault :
% * Yc - [mm], Y coordinate of the center of mass in yz plane (as function of x).
% * Iz - [mm^4], Area moment of inertia (as function of x).
% * yMtx - [mm], each line in the matrix is Y(x) from bottom to top.
% * Wxy - [mm], a matrix of the tooth width as function of Y and X.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Handle recursively a special case of TroughFaceFault %%
%%% This step considers the case of a healthy tooth, while the
%%% tooth of the other wheel is damaged.
if strcmp(defectInfo.status,'Healthy') && strcmp(defectInfo.statusOtherWheel,'ThroughFaceFault')
    defectInfo.otherWheel.status = 'Healthy' ;
    defectInfo.status = 'ThroughFaceFault' ;
    defectInfo.Y = Y ;
    invKs = calcInvShearStiff(alpha, W, G, X, Y, defectInfo) ;
    return
end % of if

%% Apply a unit force %%
F = 1 ; % unit contact force
Fs = F * cos(alpha) ; % shear component of the contact force

%% Calculate the potential energy due to shear stress Us %%
switch defectInfo.status
    case 'Healthy'
        A = W.*(2*Y) ;
        Us = calcUs(Fs, X, G, A) ;
        
    case 'ToothBreakage'
        A = defectInfo.Wx.*(2*Y) ;
        Us = calcUs(Fs, X, G, A) ;
        
    case 'PartialFaceFault'
        [Yc, Iz, A] = deal(defectInfo.geometry.Yc, defectInfo.geometry.Iz, defectInfo.geometry.A) ;
        NPtsY = defectInfo.NPtsY4Shear ;
        yMtx = calcLinspaceMtx(-Y-Yc, +Y-Yc, NPtsY) ;
        Wxy = W - defectInfo.dimensions.W*(yMtx > repmat(defectInfo.Y-Yc,1,NPtsY)) ;
        faultInds = sort([defectInfo.faultInds.start, defectInfo.faultInds.end])  ;
        faultInds = faultInds(1):faultInds(2) ;
        Us = calcUs(Fs, X, G, A, Iz, yMtx, Wxy, faultInds) ;
        
    case 'ThroughFaceFault'
        defY = defectInfo.Y ;
        A = W.*(Y+defY) ;
        Us = calcUs(Fs, X, G, A) ;
        Us = correctUForThroughFaceFault(Us, defectInfo.faultInds, Fs) ;
        
    case 'ToothDestruction'
        defY = defectInfo.Y ;
        A = W.*(Y+defY) ;
        Us = calcUs(Fs, X, G, A) ;
        
end % of switch defectInfo.status

%% Express the inverse stiffness with the potential energy %%
invKs = (2*Us)/(F^2) ;

end % of function "calcInvShearStiff"
%%
function Us = calcUs(Fs, X, G, A, Iz, yMtx, Wxy, faultInds)
%{
% Description:
% This function computes the potential strain energy resulting from
% shear stress according to beam theory.
% Improtant note: The shear force Fs may vary along X but remains
% X-independent within the cumulative integral. Consequently, the term
% involving Fs is extracted from the integral and replaced by scalar multiplication.
% =====
% Inputs & Outputs:
% * please see the description in the parent function and assist user manual.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
integAlongX = 0.6 ./(G*A) ; % initialize with the case of unifrom rectangular cross-section in the z-y plane

if nargin == 8 % inconsistent cross-section in the z-y plane
    for ind = faultInds
        Qz = cumtrapz(yMtx(ind,:), yMtx(ind,:).*Wxy(ind,:)) ;
        Qz = Qz(end) - Qz ; % boundaries of definite integral along y
        integAlongY = trapz(yMtx(ind,:), Qz.^2 ./ Wxy(ind, :)) ;
        integAlongX(ind) = integAlongY ./ (2*G*Iz(ind).^2) ;
    end % of for ind (update integAlongX)
end % of if

Us = Fs.^2.*cumtrapz(X, integAlongX) ;

end % of function 'calcUs'