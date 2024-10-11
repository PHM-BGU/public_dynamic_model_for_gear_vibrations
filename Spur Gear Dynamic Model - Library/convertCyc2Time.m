function xTime = convertCyc2Time(dCyc, xCyc, dt, rpsSig)
%{
% Description:
% This function takes any structure (vector/matrix/list) "xCyc" in the cycle
% domain (periodic as a function of the cycle), corresponding to the cycle
% vector with resolution of "dCyc", and returns the equivalent
% representation of "xCyc" (or part of it) in the time domain.
% =====
% Inputs:
% dCyc - [cyc], resolution in the cycle domain.
% xCyc - the structure x in the cycle domain
% dt - [sec], the time resolution.
% rpsSig - [rps], the rps signal (step or not) in time.
% =====
% Outputs:
% * xTime - the value of x that matches the desired cycle.
% =====
% In-Function Variables:
% * cycTime - the fractions (between 0 to 1) of the phase at time 't'
%   relatively to 2*pi.
% * cycInds - the indices in the cycle vector for each point in cycTime.
%   in case of a cycle fraction corresponding to the structure's length,
%   take the first value (due to periodicity).
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Calculate the cycle indexes %%
cycTime = cumsum(rpsSig*dt) ;
cycTime = cycTime - cycTime(1) ;
cycInds = mod(round(cycTime./dCyc), length(xCyc)) + 1 ;

%% Build xTime in time depending on the dimension of xCyc %%
if isvector(xCyc)
    xTime = xCyc(cycInds) ;
elseif length(size(xCyc)) == 2
    xTime = xCyc(:, cycInds) ;
elseif length(size(xCyc)) == 3
    xTime = xCyc(:, :, cycInds) ;
end % of if

end % of function "convertCyc2Time"