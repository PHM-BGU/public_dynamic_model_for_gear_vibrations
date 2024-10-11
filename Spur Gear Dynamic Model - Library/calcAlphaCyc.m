function alphaCyc = calcAlphaCyc(alphaNom, cycIn, tr, defectInfo)
%{
% Description:
% This function calculates the pressure angle alpha as function of the
% cycle of the input shaft.
% =====
% Inputs:
% * cycIn - [cyc/2π], cycle of the input (driving) wheel.
% * tr - transmission ratio between teeth (z_out/zIn).
% * alphaNom - [rad], nominal pressure angle.
% * defectInfo - a structure with all the information about the fault.
%   # defectInfo.status - health status of the wheel.
%   # defectInfo.faultCycs.start/end - the fault's start/end cycle point.
%   # defectInfo.otherWheel.R.base - [mm], the base raduis of the non-damaged wheel.
%   # defectInfo.defR.start/mid/end - [mm], radii from undamaged wheel center
      to the involute profile points that match the fault edges on the damaged wheel.
% =====
% Outputs:
% * alphaCyc - [rad], pressure angle as function of the input shaft cycle.
% =====
% In-Function Variables:
% * remCycIn - the remaning cycle of the in wheel after tr
% cycles, that is, after one cycle of the out wheel.
% * lenR - the length of each part of the fault (initiation-transition
%   and transition-exit).
% * [mid2startR, mid2endR] - [mm], the ranges of radii in the initiation-transition
% and transition-exit areas.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Initialize pressure angle as nominal %%
alphaCyc = alphaNom * ones(size(cycIn)) ;

%% Handle different fault cases %%
switch defectInfo.status
    case 'ThroughFaceFault'
        remCycIn = rem(cycIn, tr) ;
        [defR, Rb] = deal(defectInfo.defR, defectInfo.otherWheel.R.base) ;
        defStartInd = find(remCycIn>=defectInfo.faultCycs.start, 1, 'first') ;
        def_endInd = find(remCycIn<=defectInfo.faultCycs.end, 1, 'last') ;
        defInds = (defStartInd:def_endInd)' ;
        lenR = ceil(length(defInds)/2) ;
        mid2startR = linspace(defR.start, defR.mid, lenR) ;
        mid2endR = linspace(defR.mid, defR.end, lenR) ;
        for ii = lenR:-1:1
            alphaCyc(defInds(ii)) = alphaNom + ...
                (sqrt(defR.start^2-Rb^2)-sqrt(mid2startR(ii)^2-Rb^2))/Rb ;
            alphaCyc(defInds(ii+lenR-1))= alphaNom + ...
                (sqrt(defR.end^2-Rb^2)-sqrt(mid2endR(ii)^2-Rb^2))/Rb ;
        end % of for
end % of switch defectInfo.status
end % of function 'calcAlphaCyc'