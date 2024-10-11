function Keq = calcKeq(inWheelInfo, outWheelInfo, inWheelDefectInfo, outWheelDefectInfo)
%{
% Description:
% This function calculates the equivalent stiffness of a tooth pair
% according to the potential energy method (beam theory). The equivalent
% stiffness of a tooth pair is composed of the total stiffness of the input
% and output wheel (each includes the axial, bending, shear, and filllet
% foundation stiffness), and the hertzian contact stiffness. The
% stiffnesses are connected in series.
% =====
% Inputs:
% * (in/out)WheelInfo - a structure with all the information of the tooth.
% * (in/out)WheelDefectInfo - a structure with all the information about the fault.
% =====
% Outputs:
% * Keq - [N/mm], A vector of the equivalent stiffness of a tooth pair.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Input check - case of a healthy tooth pair %%
if nargin == 2
    [defectInfo.status, defectInfo.statusOtherWheel] = deal('Healthy') ;
    defectInfo.alpha = outWheelInfo.alpha ;
    [outWheelDefectInfo, inWheelDefectInfo] = deal(defectInfo) ;
end % of if

%% Avoid unnecessary calculations in case of a missing tooth %%
if strcmp(outWheelDefectInfo.status, 'MissingTooth') || strcmp(inWheelDefectInfo.status, 'MissingTooth')
    Keq = zeros(length(outWheelInfo.X) - outWheelInfo.initContInd + 1, 1) ;
    return
end % of if

%% Load the initial contact index and tooth width %%
inWheelDefectInfo.initContInd = inWheelInfo.initContInd ;
outWheelDefectInfo.initContInd = outWheelInfo.initContInd ;
try
    inWx = inWheelInfo.W - outWheelDefectInfo.otherWheel.defWidth ;
    outWx = outWheelDefectInfo.Wx ;
catch
    [inWx, outWx] = deal(outWheelInfo.W*ones(size(outWheelInfo.X))) ;
end % of try-catch

%% Calculate the stiffnesses composing the equivalent stiffness %%
invKIn = calcInvKWheel(inWheelInfo, inWheelDefectInfo,inWx) ;
invKOut = calcInvKWheel(outWheelInfo, outWheelDefectInfo, outWx) ;
invKh = calcInvHertzianContStiff(outWx, outWheelInfo.E, outWheelInfo.v) ;
invKh = flip(invKh(outWheelInfo.initContInd:end)) ;
Keq = 1./(invKIn + invKOut + invKh) ;

%% Handle a special case of ToothBreakage %%
%%%% similarly to the MissingTooth case, the indexes of the missing
%%%% points along the profile do not contribute stiffness to the
%%%% equivalent stiffness and therefore should be zeroed.
if strcmp(outWheelDefectInfo.status,'ToothBreakage') && outWheelDefectInfo.dimensions.tipLoss
    defHightFromTip = outWheelDefectInfo.dimensions.tipLoss ; % measurement from tooth tip [mm]
    isHealthyInd = outWheelInfo.X <= (outWheelInfo.X(end)- defHightFromTip) ;
    isHealthyInd = flip(isHealthyInd(outWheelInfo.initContInd:end)) ;
    Keq = Keq.*isHealthyInd ; % zero Keq in the indices of the missing tooth
end % of if

end % of function "calcKeq"
%%
function invKWheel = calcInvKWheel(wheelInfo, defectInfo, Wx)
%{
% Description:
% This function calculates the equivalent inverse stiffness of each wheel
% due to axial, bending, and shear stresses. In addition, the stiffness due
% to the residual stress induced by the fillet foundation is also calculated.
% =====
% Inputs:
% * wheelInfo - a structure with all the information of the tooth.
% * defectInfo - a structure with all the information about the fault.
% * Wx - [mm], the variation of the tooth width along the x axis.
% =====
% Outputs:
% * invKWheel - [N/mm]^-1, The equivalent inverse stiffness of the wheel.
% =====
% In-Function Variables:
% * invKa - The inverse axial stiffness.
% * invKs - The inverse shear stiffness.
% * invKb - The inverse bending stiffness.
% * invKf - The inverse stiffness induced by the fillet foundation.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Load all the relevant parameters %%
X = wheelInfo.X ;
Y = wheelInfo.Y ;
W = wheelInfo.W ;
E = wheelInfo.E ;
G = wheelInfo.G ;
initContInd = wheelInfo.initContInd ;
coeffMtx = wheelInfo.coeffMtx4FilletStiff ;

try
    alpha = defectInfo.alphaX ;
catch
    alpha = wheelInfo.alpha ;
end % of try-catch

%% Calculate the inverse stiffnesses according to the potential energy method %%
invKa = calcInvAxialStiff(alpha, W, E, X, Y, defectInfo) ;
invKs = calcInvShearStiff(alpha, W, G, X, Y, defectInfo) ;
invKb = calcInvBendingStiff(alpha, W, E, X, Y, defectInfo) ;

%% Calculate the inverse stiffness induced by the fillet foundation %%
invKf = calcInvFilletFoundationStiff( ...
    alpha, Wx, E, X, Y, wheelInfo.R.dedendum, wheelInfo.holeRadius, coeffMtx, defectInfo) ;

%% Calculate the equivalent inverse stiffness of the wheel %%
invKWheel = invKa + invKs + invKb + invKf ;
invKWheel = invKWheel(initContInd:end) ;

if strcmp(wheelInfo.wheel,'out')
    invKWheel = flip(invKWheel) ;
end % of if
end % of function 'calcInvKWheel'
%%