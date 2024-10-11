function plotToothFace(wheelInfo, defectInfo)
%{
% Description:
% This function plots the tooth face in the ZX plane.
% In case of a fault, the function will plot the defected tooth even if the
% other teeth are healthy.
% =====
% Inputs:
% * wheelInfo - a structure with all the wheel information (X, W, initContInd).
% * defectInfo - a structure with all the fault related information (dimensions, faultInds).
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
%}

%% Set bounds for the rectangle representing a healthy tooth face %%
[X, W] = deal(wheelInfo.X, wheelInfo.W) ;
[dx, dz] = deal(X(end)-X(1), W) ;
colorDefFace = [0.9 0.9 0.9] ;
colorTipReliefFace = [0.4 0.4 0.4] ;

%% Modify and plot tooth face According to health status %%
switch defectInfo.status
    case 'ThroughFaceFault'
        rectangle('Position',[-W/2 X(1) dz dx])
        rectangle('Position',[-W/2 X(wheelInfo.toothModification.tipRelief.edgeInd) dz (X(end)-X(wheelInfo.toothModification.tipRelief.edgeInd))], 'FaceColor', colorTipReliefFace)
        dx = X(defectInfo.faultInds.end) - X(defectInfo.faultInds.start) ;
        rectangle('Position',[-W/2 X(defectInfo.faultInds.start) dz dx], 'FaceColor', colorDefFace)
    case 'PartialFaceFault'
        rectangle('Position',[-W/2 X(1) dz dx])
        rectangle('Position',[-W/2 X(wheelInfo.toothModification.tipRelief.edgeInd) dz (X(end)-X(wheelInfo.toothModification.tipRelief.edgeInd))],'FaceColor', colorTipReliefFace)
        dx = X(defectInfo.faultInds.end) - X(defectInfo.faultInds.start) ;
        [defW, dz] = deal(defectInfo.dimensions.W) ;
        rectangle('Position',[W/2-defW X(defectInfo.faultInds.start) dz dx], 'FaceColor', colorDefFace)
    case 'MissingTooth'
        dx = X(floor(0.9*wheelInfo.initContInd)) - X(1) ;
        rectangle('Position',[-W/2 X(1) dz dx])
    case 'ToothBreakage'
        [x, z, tipLoss] = deal(defectInfo.dimensions.x, defectInfo.dimensions.z, defectInfo.dimensions.tipLoss) ;
        zPentagon = [-W/2 -W/2 (-W/2+z) +W/2 +W/2] ;
        xPentagon = [X(1)  (X(end)-x) (X(end)-tipLoss) (X(end)-tipLoss) X(1)] ;
        patch(zPentagon, xPentagon, 'w') ;
    case 'ToothDestruction'
        rectangle('Position',[-W/2 X(1) dz dx])
        rectangle('Position',[-W/2 X(wheelInfo.toothModification.tipRelief.edgeInd) dz (X(end)-X(wheelInfo.toothModification.tipRelief.edgeInd))],'FaceColor', colorTipReliefFace)
        dx = X(defectInfo.faultInds.end) - X(defectInfo.faultInds.start) ;
        rectangle('Position',[-W/2 X(defectInfo.faultInds.start) dz dx], 'FaceColor', colorDefFace)
    otherwise
        rectangle('Position',[-W/2 X(1) dz dx])
        rectangle('Position',[-W/2 X(wheelInfo.toothModification.tipRelief.edgeInd) dz (X(end)-X(wheelInfo.toothModification.tipRelief.edgeInd))],'FaceColor', colorTipReliefFace)
end % of switch-case

set(gca, 'XDir', 'reverse');
xlim([-W/2-1, W/2+1]) ;  ylim([X(1)-1, X(end)+1])
xlabel('Z [mm]') ; ylabel('X [mm]') ;

end % of fucntion 'plotToothFace'