function plotToothProfile(wheelInfo, defectInfo)
%{
% Description:
% This function plots the tooth profile in the YX plane.
% In case of a fault, the function will plot the defected tooth even if the
% other teeth are healthy.
% =====
% Inputs:
% * wheelInfo - a structure with all the wheel information (module, z, alpha, initContInd).
% * defectInfo - a structure with all the fault related information (Y, faultInds, dimensions).
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Load / Calculate tooth profile Parameters %%
[module, z, alpha, initContInd] = ...
    deal(wheelInfo.module, wheelInfo.z, wheelInfo.alpha, wheelInfo.initContInd) ;
% Notice: IF - Inactive Face ; OF - Operational Face.
[xOF, yOF, R, ~, xIF, yIF, xAddendum, yAddendum] = genHealthyToothProfile(module, z, alpha, wheelInfo.toothModification) ;
radii = [R.dedendum, R.base, R.pitch, R.addendum] ;

%% Plot characteristic circles %%
colorCharCircs = [0.75 0.75 0.75] ;
colorDefFace = [0.9 0.9 0.9] ;
centers = zeros(length(radii), 2) ;
viscircles(centers, radii, 'color', colorCharCircs, 'LineWidth', 0.4, 'LineStyle', '--') ;
hold on

%% Modify tooth profile according to health status %%
switch defectInfo.status
    case 'ThroughFaceFault'
        yOF = defectInfo.Y ;
    case 'PartialFaceFault'
        healthyY = yOF ;
        yOF = defectInfo.Y ;
        [startInd, endInd] = deal(defectInfo.faultInds.start, defectInfo.faultInds.end) ;
        fill([healthyY(startInd:endInd) ; flipud(yOF(startInd:endInd))], ...
            [xOF(startInd:endInd) ; flipud(xOF(startInd:endInd))], colorDefFace);
    case 'MissingTooth'
        ind = floor(0.9*initContInd) ;
        xOF(ind+1:end) = xOF(ind) ;
        xIF = xOF ;
        xAddendum = xOF(ind)*ones(size(xAddendum)) ;
    case 'ToothBreakage'
        [x, ~, tipLoss] = deal(defectInfo.dimensions.x, defectInfo.dimensions.z, defectInfo.dimensions.tipLoss) ;
        [~, indChipping] = min(abs(xOF-(xOF(end)-x))) ;
        [~, indLoss] = min(abs(xOF-(xOF(end)-tipLoss))) ;
        xOF(indLoss+1:end) = xOF(indLoss) ;
        xIF = xOF ;
        xAddendum = xOF(indLoss)*ones(size(xAddendum)) ;
        YChipping = [yIF(indChipping:indLoss) ; flipud(yOF(indChipping:indLoss))] ;
        XChipping = [xIF(indChipping:indLoss); flipud(xOF(indChipping:indLoss))] ;
        patch(YChipping, XChipping, XChipping, colorDefFace) ;
    case 'ToothDestruction'
        yOF = defectInfo.Y ;
        addendumInds = find(yAddendum<=yOF(end)) ;
        [xAddendum, yAddendum] = deal(xAddendum(addendumInds), yAddendum(addendumInds)) ;
end % of switch-case

%% Plot tooth profile %%
plot(yOF, xOF, yIF, xIF, yAddendum, xAddendum, 'LineWidth', 1, 'color', 'k') ;
ylim([min(xOF), max(radii)]) ;  xlim([min(yOF), max(yOF)]) ;
ylabel('X [mm]') ;  xlabel('Y [mm]') ;
axis equal

end % of function 'plotToothProfile'