function X = uniformRand(systematicErr, statisticErr, sizeX)
%{
% Description:
% This function generates a sizeX matrix of a nominal value with a desired
% random noise distributed uniformly.
% =====
% Inputs:
% * systematicErr - systematic error
% * statisticErr - statistical error
% * sizeX - size of the desired output [rows, cols]
% =====
% Outputs:
% X - the desired matrix according to the description.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by BGU-PHM Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

if length(sizeX) == 1
    sizeX = [sizeX, 1] ;
end

X = systematicErr + statisticErr.*((rand(sizeX)-0.5)*2) ;

end % of function "uniformRand"