function k = calcShaftStiff(r, L, G)
%{
% Description:
% This function calculates the torsional stiffness of a cylindrical shaft,
% based on its radius (r), length (L), and shear modulus (G).
% =====
% Inputs:
% * r [mm] - Shaft's radius.
% * L [mm] - Shaft's length.
% * G [N/mm^2] - Shear modulus of the shaft's material.
% =====
% Outputs:
% * k [Nm/rad] - Torsional stiffness of the shaft.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

k = pi * r.^4 * G ./ ( 2 * L ) ;
k = k * 1e-3 ; % from [N-mm/rad] to [Nm/rad]

end % of function 'calcShaftStiff'