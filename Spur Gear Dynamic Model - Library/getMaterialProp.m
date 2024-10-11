function [E, v, G, rho] = getMaterialProp(material)
%{
% Description:
% This function stores a basic database of relevant material properties of
% common metals used for gear manufacturing.
% =====
% Inputs:
% * material - The name of the metal.
% =====
% Outputs:
% * E - [MPa], Young modulus of elasticity.
% * v - Poisson ratio.
% * G - [MPa], Shear modulus.
% * rho - [kg/m^3], Density.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% step 1: load E,v,rho according to the desired metal %%
switch material
    case {'steel', 'Steel'}
        E = 210e3 ;
        v = 0.30 ;
        rho = 7.84e3 ;
    case {'copper', 'Copper'}
        E = 120e3 ;
        v = 0.33 ;
        rho = 8.92e3 ;
    case {'aluminum', 'Aluminum'}
        E = 70e3 ;
        v = 0.33 ;
        rho = 2.70e3 ;
    case {'brass', 'Brass'}
        E = 110e3 ;
        v = 0.36 ;
        rho = 8.70e3 ;
    case {'cast-iron', 'cast iron', 'Cast Iron', 'CastIron'}
        E = 140e3 ;
        v = 0.25 ;
        rho = 7.80e3 ;
    case {'bronze', 'Bronze'}
        E = 110e3 ;
        v = 0.35 ;
        rho = 8.30e3 ;
    otherwise
        error('The requested gear material does not exist in database.')
end % of switch-case

%% step 2: calculate shear modulus %%
G = E / 2 / (1 + v) ; 

end % of function 'getMaterialProp'