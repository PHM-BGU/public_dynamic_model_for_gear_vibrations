function KCycList = calcVariableStiffMtx(gmsCyc, cycVctr, outWheelInfo, ...
    inWheelInfo, testRigParam, outWheelDefectInfo, inWheelDefectInfo)
%{
% Description:
% This function calculates the variable stiffness matrix as a 3D array.
% The output KCycList is a list with size of ndof x ndof x lenCyc.
% The reason for calculation in cycle rather than in time is realted to
% optimization of the running time and complexity, as explained in the
% numerical solution section in detail.
% =====
% Inputs:
% * gmsCyc - [N/m] - a vector of the gms in cycle (coarse).
% * (in/out)WheelInfo - a structure with geomterical information of the tooth.
% * testRigParam -  parameters of the shafts and bearings stiffnesses.
% * (in/out)WheelDefectInfo - a structure with all the information about the fault.
% =====
% Outputs:
% * KCycList - ndof x ndof x lenCyc matrix of the structural stiffness.
% =====
% Significant In-Function Variables:
% * (in/out)Rb - [m], base radius.
% * ndof - Number degrees of freedom ( = 13).
% * defectInfo.alphaCyc - [rad], a vector of the pressure angle in cycle
%   in case of a ThroughFaceFault.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
    %}
    
    %% Extract relevant parameters, mainly geometrical %%
    [alpha, beta] = deal(outWheelInfo.alpha, outWheelInfo.beta) ;
    [outRb, inRb] = deal(outWheelInfo.R.base*1e-3, inWheelInfo.R.base*1e-3) ; % from [mm] to [m]
    [lenCyc, ndof] = deal(length(gmsCyc), 13) ;
    
    %% Calculate the variable K according to the health status %%
    switch outWheelDefectInfo.status
        case 'ThroughFaceFault'
            varKList = zeros(ndof, ndof, lenCyc) ;
            for ii = 1:lenCyc
                alphaCyc = outWheelDefectInfo.alphaCyc(ii) ;
                geomZIndep = calcGeomZIndepCoeffs(alphaCyc, beta, outRb, inRb) ;
                varKList(:, :, ii) = geomZIndep .* gmsCyc(ii) ;
            end % of for ii
        otherwise
            W = outWheelInfo.W * 1e-3 ; % from [mm] to [m]
            zMin  = -0.5*W * ones(size(gmsCyc)) ;
            zMax  = +0.5*W * ones(size(gmsCyc)) ;
            switch outWheelDefectInfo.status
                case 'ToothBreakage'
                    [~, faultStartInd] = min(abs(cycVctr - outWheelDefectInfo.defEntranceCycIn)) ;
                    [~, faultEndInd] = min(abs(cycVctr - outWheelDefectInfo.defExitCycIn)) ;
                    defWidth = outWheelDefectInfo.dimensions.z * 1e-3 ; % from [mm] to [m]
                    zMax(faultStartInd:faultEndInd) = +0.5*W - linspace(defWidth, 0, length(faultStartInd:faultEndInd))' ;
                case 'PartialFaceFault'
                    [~, faultStartInd] = min(abs(cycVctr - outWheelDefectInfo.faultCycs.start)) ;
                    [~, faultEndInd] = min(abs(cycVctr - outWheelDefectInfo.faultCycs.end)) ;
                    defWidth = outWheelDefectInfo.dimensions.W * 1e-3 ; % from [mm] to [m]
                    zMax(faultStartInd:faultEndInd) = +0.5*W - defWidth ;
            end % of switch outWheelDefectInfo.status
            geomZIndepCoeffs = calcGeomZIndepCoeffs(alpha, beta, outRb, inRb) ;
            geomZDepCoeffs = calcGeomZDepCoeffs(alpha, beta, outRb, inRb, zMin, zMax)  ;
            geomCoeffs = geomZIndepCoeffs + geomZDepCoeffs ;
            geomCoeffs = reshape(geomCoeffs,[ndof^2, lenCyc]) ;
            varKVctr = geomCoeffs.*repmat(gmsCyc.', ndof^2, 1) ;
            varKList = reshape(varKVctr,[ndof, ndof, lenCyc]);
    end % of switch outWheelDefectInfo.status
    
    %% Consider z-dependency in the geometry matrices %%
    constK = calcConstStiffMtx(testRigParam, inWheelInfo.rotStiff, outWheelInfo.rotStiff) ;
    KCycList = varKList + constK ;
    
end % of function 'calcVariableStiffMtx'
%%
function constK = calcConstStiffMtx(testRigParam, inRotStiff, outRotStiff)
%{
% Description:
% This function calculates the constant diagonal components of the
% stiffness matrix, according to the following coordinates vector:
% u = [xOut, yOut, xIn, yIn, thetaOut, thetaIn, thetaBrake, zOut, zIn, phiOut, phiIn, psiOut, psiIn]
% The only coupling in the constant K is that of the output wheel
% and the brake, expressed by the torsional stiffness of the output shaft.
% =====
% Inputs:
% * testRigParam - parameters of the shafts and bearings stiffnesses:
    # testRigParam.bearingStiff - [N/m], bearing stiffness.
    # testRigParam.shaftPolarStiff.(wheel) - [Nm] shaft torsion stiffness.
    # (in/out)rotStiff - [Nm] - rotational stiffness of the rigid
      cylinders (that is, the wheels).
% =====
% Outputs:
% * constK - The constant diagonal stiffness matrix.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
bearingStiff = testRigParam.bearingStiff ;
inShaftPolarStiff = testRigParam.shaftPolarStiff.in ;
outShaftPolarStiff = testRigParam.shaftPolarStiff.out ;

constK = diag([bearingStiff, bearingStiff, bearingStiff, ...
    bearingStiff, outShaftPolarStiff, inShaftPolarStiff, ...
    outShaftPolarStiff, bearingStiff, bearingStiff, ...
    outRotStiff, inRotStiff, outRotStiff, inRotStiff]) ;

[constK(5,7), constK(7,5)] =  deal(- outShaftPolarStiff) ;
end % of function 'calcConstStiffMtx'
%%
function geomZDepCoeffs = calcGeomZDepCoeffs(alpha, beta, outRb, inRb, zMin, zMax)
%{
% Description:
% This function calculates a matrix of geometrical coefficients of
% components in the stiffness matrix that are z-dependent.
% Since the gms is assumed to distribute uniformly anlong the z-axis,
% the mean coordinate of the tooth width and the squared z-axis are
% calculated and considered seperately.
% This separation is meant to improve efficiency and running time.
% More information is detailed in the manual.
% =====
% Inputs:
% * alpha - [rad], pressure line (nominal or instantaneous).
% * beta - [rad], helix angle (zero for spur gears).
% * (in/out)Rb - [m], base radius.
% * [zMin, zMax] - [m], vectors with the edges of the tooth width.
%   The length of the vector is the length of the cycle vector.
% =====
% Outputs:
% * geomZDepCoeffs - an ndof x ndof x lenCyc matrix.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
%% Call coefficients only relevant for the z-dependent geometry %%
[~, ~, ~, ~, ~, ~, ~, f8, ~, ~, ~, ~, f13, ~, ~, ~, ~, f18, f19, ...
    ~, ~, ~, f23, f24] = calcFCoeffs(alpha, beta, outRb, inRb) ;

%% Build the coefficient matrix: [topLeft, topRight ; botLeft, botRight] %%
topLeft = zeros(9) ;
topRight =  [...
    +f8.*f23, -f8.*f23, +f8.*f24, -f8.*f24 ; ...
    +f13.*f23, -f13.*f23, +f13.*f24, -f13.*f24 ; ...
    -f8.*f23, +f8.*f23, -f8.*f24, +f8.*f24 ; ...
    -f13.*f23, +f13.*f23, -f13.*f24, +f13.*f24 ; ...
    +f18.*f23, -f18.*f23, +f18.*f24, -f18.*f24 ; ...
    +f19.*f23, -f19.*f23, +f19.*f24, -f19.*f24 ; ...
    0, 0, 0, 0 ; ...
    +sin(beta).*f23, -sin(beta).*f23, +sin(beta).*f24, -sin(beta).*f24 ; ...
    -sin(beta).*f23, +sin(beta).*f23, -sin(beta).*f24, +sin(beta).*f24 ...
    ] ;
botLeft = topRight' ;
botRight =  [...
    +f23.*f23, -f23.*f23, +f23.*f24, -f23.*f24 ; ...
    -f23.*f23, +f23.*f23, -f24.*f23, +f23.*f24 ; ...
    +f23.*f24, -f24.*f23, +f24.*f24, -f24.*f24 ; ...
    -f23.*f24, +f23.*f24, -f24.*f24, +f24.*f24 ; ...
    ] ;
geomZDepCoeffs = [topLeft, topRight ; botLeft, botRight] ;

%% Create a list varying in cycle %%
lenCyc = length(zMin) ;
zMean = reshape( (zMin + zMax)/2, [1, 1, lenCyc] ) ;
zSquaredMean = reshape( (zMin.^2 + zMin.*zMax + zMax.^2)/3, [1,1,lenCyc] ) ;
geomZDepCoeffs = repmat( geomZDepCoeffs, 1, 1, lenCyc) ;
expectationMtx = [zeros(9,9,lenCyc), zMean.*ones(9,4) ; zMean.*ones(4,9), zSquaredMean.*ones(4,4,lenCyc)] ;
geomZDepCoeffs = geomZDepCoeffs .* expectationMtx ;

end % of function 'calcGeomZDepCoeffs'
%%
function geomZIndepCoeffs = calcGeomZIndepCoeffs(alpha, beta, outRb, inRb)
%{
% Description:
% This function calculates a matrix of geometrical coefficients of
% components in the stiffness matrix that are z-independent.
% This separation is meant to improve efficiency and running time.
% More information is detailed in the manual.
% =====
% Inputs:
% * alpha - [rad], pressure line (nominal or instantaneous).
% * beta - [rad], helix angle (zero for spur gears).
% * (in/out)Rb - [m], base radius.
% =====
% Outputs:
% * geomZIndepCoeffs - an ndof x ndof matrix.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}
%% Call coefficients only relevant for the z-independent geometry %%
[f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, ...
    f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f22] = ...
    calcFCoeffs(alpha, beta, outRb, inRb) ;

%% Construct only the upper triangle of the coefficient matrix (due to symmetry) %%
geomUpperTri = [...
    +f5, f6, -f5, -f6, +f9, +f10, 0, +f7, -f7, +f8.*f1, +f8.*f3, +f8.*f2, +f8.*f4 ; ...
    0, +f11, -f6, -f11, +f14, +f15, 0, +f12, -f12, +f13.*f1, +f13.*f3, +f13.*f2, +f13.*f4 ; ...
    0, 0, +f5, +f6, -f9, -f10, 0, -f7, +f7, -f8.*f1, -f8.*f3, -f8.*f2, -f8.*f4 ; ...
    0, 0, 0, +f11, -f14, -f15, 0, -f12, +f12, -f13.*f1, -f13.*f3, -f13.*f2, -f13.*f4 ; ...
    0, 0, 0, 0, +f20, +f21, 0, +f16, -f16, +f18.*f1, +f18.*f3, +f18.*f2, +f18.*f4 ; ...
    0, 0, 0, 0, 0, +f22, 0, +f17, -f17, +f19.*f1, +f19.*f3, +f19.*f2, +f19.*f4 ; ...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; ...
    0, 0, 0, 0, 0, 0, 0, +sin(beta).^2, -sin(beta).^2, +sin(beta).*f1, +sin(beta).*f3, +sin(beta).*f2, +sin(beta).*f4 ; ...
    0, 0, 0, 0, 0, 0, 0, 0, +sin(beta).^2, -sin(beta).*f1, -sin(beta).*f3, -sin(beta).*f2, -sin(beta).*f4 ; ...
    0, 0, 0, 0, 0, 0, 0, 0, 0, +f1.*f1, +f1.*f3, +f1.*f2, +f1.*f4 ; ...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, +f3.*f3, +f2.*f3, +f3.*f4 ; ...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, +f2.*f2, +f2.*f4 ; ...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, +f4.*f4 ...
    ] ;

%% Build the coefficient matrix by manipulations on the upper triangle %%
geomZIndepCoeffs = geomUpperTri + geomUpperTri' - diag(diag(geomUpperTri)) ;

end % of function 'calcGeomZIndepCoeffs'
%%
function [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, ...
    f15, f16, f17, f18, f19, f20, f21, f22, f23, f24] = ...
    calcFCoeffs(alpha, beta, outRb, inRb)
%{
% Description:
% This function calculates 24 geometrical terms that are later used for the
calculation of the geometrical coefficient matrix of the stiffness matrix.
% =====
% Inputs:
% * alpha - [rad], pressure line (nominal or instantaneous).
% * beta - [rad], helix angle (zero for spur gears).
% * (in/out)Rb - [m], base radius.
% =====
% Outputs:
% * fi (i=1,2,...,24) - geometrical terms.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
    %}
    %% Calculate the pressure line length %%
    [outL, inL] = deal(outRb*tan(alpha), inRb*tan(alpha)) ;
    
    %% Calculate geometrical terms prior to matrices formulation %%
    f1 = -outRb.*sin(alpha).*sin(beta) + cos(alpha).*outL.*sin(beta) ;
    f2 = -(outRb.*cos(alpha) + outL.*sin(alpha)).*sin(beta) ;
    f3 = -inRb.*sin(alpha).*sin(beta) + cos(alpha).*inL*sin(beta) ;
    f4 = -inRb.*cos(alpha).*sin(beta) - sin(alpha).*inL*sin(beta) ;
    f5 = cos(beta).^2.*sin(alpha).^2 ;
    f6 = cos(alpha).*cos(beta).^2.*sin(alpha) ;
    f7 = cos(beta).*sin(alpha).*sin(beta) ;
    f8 = cos(beta).*sin(alpha) ;
    f9 = outRb.*cos(beta).^2.*sin(alpha) ;
    f10 = inRb.*cos(beta).^2.*sin(alpha) ;
    f11 = cos(alpha).^2.*cos(beta).^2 ;
    f12 = cos(alpha).*cos(beta).*sin(beta) ;
    f13 = cos(alpha).*cos(beta) ;
    f14 = outRb.*cos(alpha).*cos(beta).^2 ;
    f15 = inRb.*cos(alpha).*cos(beta).^2 ;
    f16 = outRb.*cos(beta).*sin(beta) ;
    f17 = inRb.*cos(beta).*sin(beta) ;
    f18 = outRb.*cos(beta) ;
    f19 = inRb.*cos(beta) ;
    f20 = outRb.^2.*cos(beta).^2 ;
    f21 = outRb.*inRb.*cos(beta).^2 ;
    f22 = inRb.^2.*cos(beta).^2 ;
    f23 = -cos(alpha) ;
    f24 = +sin(alpha) ;
end % of function 'calcFCoeffs'