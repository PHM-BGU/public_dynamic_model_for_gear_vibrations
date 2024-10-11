function xMtx = calcLinspaceMtx(starts, ends, NPts)
%{
% Description:
% Creates a matrix xMtx, where each row is a vector with linear spacing
% between the corresponding values in the starts and ends vectors.
% =====
% Inputs:
% * starts - Vector of the start values.
% * ends - Vector of the end values.
% * NPts - Number of points to generate for each range.
% =====
% Outputs:
% * xMtx - Matrix where each row is a linspace vector.
% =====
% In-Function Variables:
% * deltaVctr - the vector of resolutions
% * deltaMtx - vector of steps
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Input check %%
if ~iscolumn(starts)
    starts = starts' ;
end % of if

if ~iscolumn(ends)
    ends = ends' ;
end % of if

%% Build the linear space matrix X %%
deltaVctr = (ends-starts) / (NPts-1) ; % resolutions vector
deltaMtx = deltaVctr * [0:(NPts-1)] ; % matrices multiplication
xMtx = repmat(starts, 1, NPts) + deltaMtx ;

end % of function 'calcLinspaceMtx'