function invKb = calcInvBendingStiff(alpha, W, E, X, Y, defectInfo)
%{
% Description:
% This function calculates the inverse stiffness derived from the
% potential energy of the tooth due to bending stresses.
% The energy is calculated by integration using the cumtrapz function.
% The function performs non-trivial calculations that are based on
% mathematical and physical analysis, described thoroughly in the manual.
% =====
% Inputs:
% * alpha - [rad], pressure angle, may be a vector.
% * W - [mm], tooth width.
% * E - [MPa], Young modulus of Elasticity.
% * [X,Y] - [mm], Cartesian coordinates of a healthy tooth profile.
% * defectInfo - a structure with relevant information about the fault.
% =====
% Outputs:
% * invKb - [N/mm]^-1, A vector of the inverse bending stiffness.
% =====
% In-Function Variables:
% * Ub - potential energy due to bending stresses.
% * [Fs,Fa] - shear and axial components of the contact force, respectively
% * [Yc, Zc] - [mm], Y & Z coordinate of the center of mass in yz plane (as function of x).
% * [Iz, Iy, Iyz] - [mm^4], Area moments of inertia (as function of x).
% * ZRsltntFa - [mm], The radius from the origin to the resultant axial force (Fa) in xz plane.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Handle recursively cases of different faults %%
%%% This step considers the case of a healthy tooth, while the
%%% tooth of the other wheel is damaged.
if strcmp(defectInfo.status,'Healthy') && ~strcmp(defectInfo.statusOtherWheel,'Healthy')
    switch defectInfo.statusOtherWheel
        case 'ToothBreakage'
            defectInfo.statusOtherWheel = 'Healthy' ;
            defectInfo.status = 'ToothBreakage' ;
            invKb = calcInvBendingStiff(alpha, W, E, X, Y, defectInfo) ;
        case 'PartialFaceFault'
            defectInfo.statusOtherWheel = 'Healthy' ;
            faultInds = sort([defectInfo.faultInds.start, defectInfo.faultInds.end]) ;
            faultInds = faultInds(1) : faultInds(2) ;
            WCont = W*ones(size(X)) ; % initialize
            WCont(faultInds) = WCont(faultInds) - defectInfo.otherWheel.dimensions.W ;
            invKb = calcInvBendingStiff(alpha, WCont, E, X, Y, defectInfo) ;
        case 'ThroughFaceFault'
            defectInfo.statusOtherWheel = 'Healthy' ;
            defectInfo.status = 'ThroughFaceFault' ;
            defectInfo.Y = Y ;
            invKb = calcInvBendingStiff(alpha, W, E, X, Y, defectInfo) ;
        otherwise
            defectInfo.statusOtherWheel = 'Healthy' ;
            invKb = calcInvBendingStiff(alpha, W, E, X, Y, defectInfo) ;
    end % of switch 'defectInfo.statusOtherWheel'
    return
end % of if

%% Apply a unit force %%
F = 1 ; % unit contact force
[Fs, Fa] = deal(F*cos(alpha), F*sin(alpha)) ; % shear & axial components of the contact force

%% Calculate the potential energy due to bending stress Ub %%
switch defectInfo.status
    case 'Healthy'
        Iz = W.*((2*Y).^3) / 12 ;
        Ub = calcUb(X, Y, E, Fs, Fa, Iz) ;
        
    case 'ToothBreakage'
        Wx = defectInfo.Wx ;
        Iz = Wx.*((2*Y).^3) / 12 ;
        Iy = (2*Y).*Wx.^3 / 12 ;
        [Iyz, Yc] = deal(0) ;
        ZRsltntFa = 0.5*defectInfo.defWidth ;
        Ub = calcUb(X, Y, E, Fs, Fa, Iz, Iy, Iyz, Yc, defectInfo.Zc, ZRsltntFa) ;
        
    case 'PartialFaceFault'
        faultInds = defectInfo.faultInds.start : defectInfo.faultInds.end ;
        defW = zeros(size(X)) ; % initialize
        defW(faultInds) = defectInfo.dimensions.W ;
        [Zc, Yc, Iz, Iy, Iyz] = deal(defectInfo.geometry.Zc, ...
            defectInfo.geometry.Yc, defectInfo.geometry.Iz, ...
            defectInfo.geometry.Iy, defectInfo.geometry.Iyz) ;
        ZRsltntFa = 0.5*defW ;
        Ub = calcUb(X, Y, E, Fs, Fa, Iz, Iy, Iyz, Yc, Zc, ZRsltntFa) ;
        
    case 'ThroughFaceFault'
        Iz = W.*((Y + defectInfo.Y).^3) / 12 ;
        Iy = 1 ; % assigning a neutral value to save calculations.
        Yc = 0.5*(defectInfo.Y - Y) ;
        [ZRsltntFa, Zc, Iyz]  = deal(0) ;
        Ub = calcUb(X, Y, E, Fs, Fa, Iz, Iy, Iyz, Yc, Zc, ZRsltntFa, defectInfo.faultInds) ;
        
    case 'ToothDestruction'
        Iz = W.*((Y + defectInfo.Y).^3) / 12 ;
        Iy = 1 ; % assigning a neutral value to save calculations.
        Yc = 0.5*(defectInfo.Y - Y) ;
        [ZRsltntFa, Zc, Iyz]  = deal(0) ;
        Ub = calcUb(X, Y, E, Fs, Fa, Iz, Iy, Iyz, Yc, Zc, ZRsltntFa) ;
        
end % switch defectInfo.status

%% step 4: express the inverse stiffness with the potential energy %%
invKb = (2*Ub)/(F^2) ;

end % of function "calcInvBendingStiff"
%%
function Ub = calcUb(X, Y, E, Fs, Fa, Iz, Iy, Iyz, Yc, Zc, ZRsltntFa, faultInds)
%{
% Description:
% This function is designed to compute the potential energy resulting from
% bending stress. The process involves dividing the calculations into three
% separate terms, which helps to optimize the running time, as outlined in
% the user manual. It should be noted that if the neutral axis intersects
% the origin, then the values for Iyz and My can be assumed to be zero.
% In this scenario, the calculations are significantly simplified,
% making them easier to compute.
% This function computes the potential strain energy resulting from
% bending stresses along y and z axes according to beam theory.
% Improtant note: The axial (Fa) and shear (Fs) forces, and the distance
% ZRsltntFa may vary along X but remains X-independent within the cumulative
% integral. Consequently, the terms involving them are extracted from
% the integral and replaced by scalar multiplication.
% =====
% Inputs & Outputs:
% *  please see the description in the parent function.
% =====
% In-Function Variables:
% * D_$$ - [N-m^2], flexural_rigidity. A property from beam theory that
% measures a material's resistance to bending deformation in different
% modes (y, z, yz). Calculated for each mode in the energy formula.
% * Ub_$$_$, Ub_F$F$ - Please assist user manual for explanasions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Input check (relevant for the healthy case) %%
if nargin == 6
    Iy = 1 ; % assigning a neutral value to save unnecessary calculations.
    [Iyz, Yc, Zc, ZRsltntFa] = deal(0) ;
end % of if

[Ub_FaFs, Ub_FaFa, Ub_FsFs] = deal(zeros(size(X))) ;

%% Calculate the strain energy in mixed-mode beding (yz) %%
% Clarification: mixed-mode beding (yz) applies only if Iyz is not zero %
if any(Iyz)
    Dyz = E*(Iy.*Iz - Iyz.^2)./Iyz ;
    Dz = E*(Iy.*Iz - Iyz.^2)./Iy ;
    Dy = E*(Iy.*Iz - Iyz.^2)./Iz ;
    
    Ub_yz_1 = Fa.*Fs.*ZRsltntFa.*cumtrapz(X, +X./Dyz) ;
    Ub_yz_2 = Fa.*Fs.*ZRsltntFa.*X.*cumtrapz(X, -1./Dyz) ;
    Ub_yz_3 = Fa.^2.*ZRsltntFa.*Y.*cumtrapz(X, +1./Dyz) ;
    Ub_yz_4 = Fa.^2.* ZRsltntFa.*cumtrapz(X, -Yc./Dyz) ;
    Ub_yz_5 = Fa.*Fs.*cumtrapz(X, -Zc.*X./Dyz) ;
    Ub_yz_6 = Fa.*Fs.*X.*cumtrapz(X, +Zc./Dyz) ;
    Ub_yz_7 = Fa.^2.*Y.*cumtrapz(X, -Zc./Dyz) ;
    Ub_yz_8 = Fa.^2.*cumtrapz(X, Zc.*Yc./Dyz) ;
    
    Ub_FaFs = Ub_FaFs + Ub_yz_1 + Ub_yz_2 + Ub_yz_5 + Ub_yz_6 ;
    Ub_FaFa = Ub_FaFa + Ub_yz_3 + Ub_yz_4 + Ub_yz_7 + Ub_yz_8 ;
else
    Dz = E*Iz ; % simpified term for this case
    Dy = E*Iy ;
end % of if

%% Calculate the strain energy derived from bending along z-axis %%
% Clarification: This mode applies always and for all health statuses %
Ub_z_1 =  Fs.^2.*cumtrapz(X, +X.^2./Dz/2) ;
Ub_z_2 =  Fa.*Fs.*cumtrapz(X, -X.*Yc./Dz) ;
Ub_z_3 =  Fa.^2.*cumtrapz(X, +Yc.^2./Dz/2) ;
Ub_z_4 =  Fs.^2.*X.*cumtrapz(X, -X./Dz) ;
Ub_z_5 =  Fa.*Fs.*Y.*cumtrapz(X, +X./Dz) ;
Ub_z_6 =  Fa.*Fs.*X.*cumtrapz(X, +Yc./Dz) ;
Ub_z_7 =  Fa.^2.*Y.*cumtrapz(X, -Yc./Dz) ;
Ub_z_8 =  Fs.^2.*X.^2.*cumtrapz(X, +1./Dz/2) ;
Ub_z_9 =  Fa.*Fs.*X.*Y.*cumtrapz(X, -1./Dz) ;
Ub_z_10 = Fa.^2.*Y.^2.*cumtrapz(X, +1./Dz/2) ;

Ub_FsFs = Ub_FsFs + Ub_z_1 + Ub_z_4 + Ub_z_8 ;
Ub_FaFs = Ub_FaFs + Ub_z_2 + Ub_z_5 + Ub_z_6 + Ub_z_9 ;
Ub_FaFa = Ub_FaFa + Ub_z_3 + Ub_z_7 + Ub_z_10 ;

%% Calculate the strain energy derived from bending along y-axis %%
% Clarification: This mode applies only if z-axis symmetry is distorted %
if any(ZRsltntFa) || any(Zc)
    Ub_y_1 = Fa.^2.*ZRsltntFa.^2.*cumtrapz(X, 1./Dy/2) ;
    Ub_y_2 = Fa.^2.*ZRsltntFa.*cumtrapz(X, -Zc./Dy) ;
    Ub_y_3 = Fa.^2.*cumtrapz(X, Zc.^2./Dy/2) ;
    
    Ub_FaFa = Ub_FaFa + Ub_y_1 + Ub_y_2 + Ub_y_3 ;
end % of if

%% Sum-up the strain energies from all bending modes %%
if length(unique(Fa))>1 % meaning, pressure angle is variable in X
    Ub_FaFa = correctUForThroughFaceFault(Ub_FaFa, faultInds, Fa) ;
    Ub_FaFs = correctUForThroughFaceFault(Ub_FaFs, faultInds, Fa, Fs) ;
    Ub_FsFs = correctUForThroughFaceFault(Ub_FsFs, faultInds, Fs) ;
end % of if

Ub = Ub_FaFa + Ub_FaFs + Ub_FsFs ;

end % of function 'calcUb'