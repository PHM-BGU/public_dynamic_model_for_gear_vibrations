function invKh = calcInvHertzianContStiff(Wx,E,v)
%{
% Description:
% This function calculates the inverse stiffness derived from the elastic
% force and deflection induced by the Hertzian contact stress.
% =====
% Inputs:
% Wx - [mm], tooth width - vector along the tooth axis (x axis).
% E - [MPa], Young's Modulus.
% v - Poisson's ratio.
% =====
% Outputs:
% invKh - [N/mm]^-1, The inverse of the Hertzian contact stiffness.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
Kh = pi * E * Wx / (4 * (1 - v^2 )) ;
invKh = 1./Kh ;
end % of function "calcInvHertzianContStiff"