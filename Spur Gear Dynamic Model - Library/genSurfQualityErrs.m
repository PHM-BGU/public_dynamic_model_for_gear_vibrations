function wheelInfo = genSurfQualityErrs(wheelInfo)
%{
% Description:
% This function takes the parameters of the simulated surface roughness 
% as input and generates corresponding errors in the teeth profile, unless
% there is already a pre-defined error matrix. The profile error of each 
% tooth is generated seperately and independently. 
% The profile error is the deviation from the ideal involute profile,
% meaning that the error of an ideal profile is 0.
% *Important note*: This function supports only DIN-3962 standard, that is,
% Tolerances for Cylindrical Gear Teeth (1980).
% =====
% Inputs:
% * wheelInfo - the relevant information about the wheels:
%   # wheelInfo.surfQuality - [int], tooth surface quality grade (DIN-3962).
%   # wheelInfo.z - number of teeth.
%   # wheelInfo.R.base - [mm], base radius.
%   # wheelInfo.(X/Y) - [mm], tooth profile coordinates.
%   # wheelInfo.surfQualityErrs.NPts - coarse resolution of the profile errors.
%   # wheelInfo.surfQualityErrs.statErrTyp - [%], typical statistical 
%    error of all the deviation parameters.
% =====
% Outputs:
% * wheelInfo.surfQualityErrs.errsMtx - [mm], A matrix with the profile
%   errors vector for each tooth in each column.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Return in case of an ideal profile or if an error matrix already exists %%
if isfield(wheelInfo.surfQualityErrs, 'errsMtx') % existing error matrix
    return
elseif ~wheelInfo.surfQuality % ideal profile
    wheelInfo.surfQualityErrs.errsMtx = ...
        zeros(wheelInfo.surfQualityErrs.NPts, wheelInfo.z) ;
    return
end % of if

%% Get all the relevant geometric parameters %%
surfQuality = wheelInfo.surfQuality ;
initContInd = wheelInfo.initContInd ;
baseR = wheelInfo.R.base ;
z = wheelInfo.z ;
NPts = wheelInfo.surfQualityErrs.NPts ;
initPt = norm([wheelInfo.X(initContInd), wheelInfo.Y(initContInd)]) ;
finalPt = norm([wheelInfo.X(end), wheelInfo.Y(end)]) ;
RContPts = linspace(initPt, finalPt, NPts).' ;
sContPts =  baseR * tan(acos(baseR./RContPts)) ; % the distance s = sqrt(Rp^2-Rb^2)

%% Get all the parameters of the flank deviations formula %%
[f_Ha, f_fa] = getDIN3962Tolerances(surfQuality) ;
[noiseCoeff, nCyc] = getComplementaryParameters(surfQuality, z) ;
statErr = wheelInfo.surfQualityErrs.statErrTyp ; % in range of [0,1]
f_Ha = buildDeviationsMtx(f_Ha, statErr*f_Ha, z, NPts) ;
f_fa = buildDeviationsMtx(f_fa, statErr*f_fa, z, NPts) ;
nCyc = buildDeviationsMtx(nCyc, statErr*nCyc, z, NPts) ;
phi = buildDeviationsMtx(0, pi, z, NPts) ;
sHat = repmat(sContPts-sContPts(1), 1, z) / (sContPts(end)-sContPts(1)) ;
noise = uniformRand(0, noiseCoeff, [NPts, z]) ;

%% Calculate flank deviations according to formula %%
surfQualityErrsMtx = ...
    f_Ha.*sHat + ...
    0.5*f_fa.*sin(2*pi*nCyc.*sHat + phi) + ...
    noise ;

%% Consider the wheel's role (driving/driven) %%
%{
The profile error vector of each tooth begins at the bottom of the 
involute, near the base, and extends to the tooth tip.
In this model, a new engagement starts at the bottom of the driving wheel 
and at the tip of the driven wheel.
Consequently, the errors of the driven output wheel are flipped to align 
with the contact along the pressure line
______________________________
%%%%%%___%%\          /%%%%%%%%%
%%%%%/    \%\%%GEAR%%/%%%%%%%%%%
%%%%/      \%\      /%%%%%%%%%%%
%%%/%PINION%\O\____/%%%%%%%%%%%%
%%/          \%%%%%%%%%%%%%%%%%%
________________________________
%}
if strcmp(wheelInfo.wheel, 'out')
    surfQualityErrsMtx = flipud(surfQualityErrsMtx) ;
end % of if

%% Convert from [μm] to [mm]
wheelInfo.surfQualityErrs.errsMtx = surfQualityErrsMtx * 1e-3 ;
wheelInfo.surfQualityErrs.errsMtxUnits = '[mm]' ;

end % of function 'genSurfQualityErrs'
%%
function [f_Ha, f_fa] = getDIN3962Tolerances(surfQuality)
%{
% Description:
% This function holds data of flank deviations according to DIN-3962
% DIN-3962 standard describes the Tolerances for Cylindrical Gear Teeth (1980)
% =====
% Inputs:
% * surfQuality - [int], precision grade. for example, an
% input of surfQuality=7 corresponds to DIN-7 quality.
% =====
% Outputs:
% f_Ha - [μm], profile angle deviation
% f_fa - [μm], profile form deviation
%}
f_HaVctr = [1, 1.5, 2, 3, 4.5, 6, 9, 12, 18, 28, 45, 71] ; % [μm]
f_faVctr = [1.5, 2, 3, 4, 6, 8, 11, 16, 22, 36, 56, 90] ; % [μm]
f_Ha = f_HaVctr(surfQuality) ;
f_fa = f_faVctr(surfQuality) ;
end % of function 'getDIN3962Tolerances'
%%
function [noiseCoeff, nCyc] = getComplementaryParameters(surfQuality, z)
%{
% Description:
% This function holds estimated values of the deterministic and random
% components in the formula of the manufacturing profile errors.
% =====
% Inputs:
% * surfQuality - [int], precision grade.
% * z - [int], number of teeth.
% =====
% Outputs:
% noiseCoeff - [um], coefficient of the random noise component.
% nCyc - number of cycles of the sine function
    %}
    noiseCoeffs = [2, 4, 6, 8, 10, 20, 30, 40, 60, 100, 150, 230] / 100 ;
    noiseCoeff = noiseCoeffs(surfQuality) ;
    nCyc = 2.5 + 1*(z>=20) ; % set a rule of thumb
end % of function 'getComplementaryParameters'
%%
function mtx = buildDeviationsMtx(nomValue, statErr, z, N)
%{
% Description:
% This function takes a desired nominal value and returns a N×z matrix of
% the nominal value, with random noise distributed uniformly with statErr.
% =====
% Inputs:
% * nomValue - nominal value (a scalar).
% * statErr - statistical error around the nominal value.
% * z - number of coloums in the matrix
% * N - number of rows in the matrix
% =====
% Outputs:
% mtx - N×z matrix according to the description.
%}
vctr = uniformRand(nomValue, statErr, z).' ;
mtx  = repmat(vctr, N, 1) ;
end % of function 'buildDeviationsMtx'