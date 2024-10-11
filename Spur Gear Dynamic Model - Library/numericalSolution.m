function [u, u_dot, u_ddot] = numericalSolution( ...
    EulerLagrangeComponents, initialConds, simParam)
%{
% Description:
% This function solves numerically the time-variable set of differential
% equations of motion from Euler-Lagrange model according to Newmark method
% with additional stage of Newton-Raphson method. The numerical solution is
% designed to accelerate convergence and reduce running time, as explained
% thoroughly in the user manual.
% =====
% Inputs:
% * EulerLagrangeComponents - a structure with the follwing fields:
%   # M - generalized diagonal mass matrix [ndof x ndof]
%   # KCycList - variable stiffness matrix in cycle [ndof x ndof x cycle].
%   # C - modal damping matrix (assuming uniform damping ratio).
%   # FexVctrT - excitation forces vector in time [ndof x t] .
% * initialConds - vectors with the initial conditions:
%   # initialConds.displacement/velocity - u(t=0) / u_dot(t=0).
% * simParam - a structure with relevant information for the numerical
% solution and the convergence criterion:
%   # simParam.Fs - [Hz] - sampling rate.
%   # simParam.dCycGMSCoarse - cycle resolution that is used to find
%     the corresponding stiffness matrix in each time step.
%   # simParam.numSol.beta/gamma - parameters for Newmark method.
%   # simParam.numSol.convergenceCriterion - for NR method [0-1].
%   # simParam.numSol.maxIterationsNR - maximal number of iterations of NR.
% =====
% Outputs:
% * u - matrix of displacement in time for all coordinates.
% * u_dot - matrix of velocity in time for all coordinates.
% * u_ddot - matrix of acceleration in time for all coordinates.
% =====
% Significant In-Function Variables:
% * cycIn - the rotational cycle of the driving wheel.
% * [a,b] - parameters of Newmark method, see manual for further details.
% * KCycModified - K matrix according to Newmark method.
% * invKCycModified - the inverse K matrix, calculated seperately in
% order to improve running time.
% * deltaFModified - the update of the force according to Newmark method.
% * P - the update for u, relaxed by Newton-Raphson method.
% * fP - the force according to Newton-Raphson method.
% * dP - the update in the current Newton-Raphson iteration.
% * [du, du_dot, du_ddot] - update of the current time step.
% =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by PHM-BGU Laboratory, Ben-Gurion University of the Negev,
% Be'er Sheva, Israel. 2023.
    %}
    
    %% Load Euler-Lagrange components %%
    [M, C, KCyc, FexVctrT] = deal(...
        EulerLagrangeComponents.M,...
        EulerLagrangeComponents.C,...
        EulerLagrangeComponents.KCycList,...
        EulerLagrangeComponents.FexVctrT) ;
    [Fs, dCycGMS] = deal(simParam.Fs, simParam.dCycGMSCoarse) ;
    dt = 1/Fs ;
    [ndof, lenT] = size(FexVctrT) ;
    lenCyc = length(KCyc) ;
    
    %% Pre-allocation & initial conditions %%
    [u_ddot, u_dot, u] = deal(zeros(size(FexVctrT))) ;
    [u_dot(:,1), u(:,1)] = deal(initialConds.velocity, initialConds.displacement) ;
    
    %% Initial calculations for Newmark's method %%
    [beta, gamma] = deal(simParam.numSol.beta, simParam.numSol.gamma) ;
    a = (1/(beta*dt))*M + (gamma/beta)*C ;
    b = (1/(2*beta))*M + dt*((gamma/(2*beta))-1)*C ;
    
    cycIn = u(6, 1) / (2*pi) ; % u(6,:) is the rotational angle of the driving wheel.
    cycInd = mod(round(cycIn/dCycGMS), lenCyc) + 1 ;
    u_ddot(:,1) = inv(M)*(FexVctrT(:,1) - C*u_dot(:,1)- KCyc(:, :, cycInd)*u(:,1)) ;
    
    %% Prepare inverse modified K matrix for computetional efficiency %%
    KCycModified = KCyc + repmat(((gamma/(beta*dt))*C + (1/(beta*dt^2))*M), 1, 1, lenCyc) ;
    P0 = zeros(ndof,1) ;
    invKCycModified = zeros(size(KCycModified)) ;
    for ii = 1:lenCyc
        invKCycModified(:, :, ii) = inv(KCycModified(:, :, ii)) ;
    end % of for ii
    
    %% Solve equations iteratively over time %%
    if simParam.numSol.waitBar
        h = waitbar(0, 'Processing...') ;
    end % of if
    
    for ii = 1:lenT-1
        deltaFex = FexVctrT(:, ii+1) - FexVctrT(:, ii) ;
        deltaFModified = deltaFex + a*u_dot(:, ii) + b*u_ddot(:, ii) ;
        
        %% Perform Newton-Raphson to improve the estimation of P %%
        cycIn = u(6, ii) / (2*pi) ;
        cycInd = mod(round(cycIn/dCycGMS), lenCyc) + 1 ;
        KModified  = KCycModified(:, :, cycInd) ;
        P = P0 ;
        KModified_0 = KModified ;
        fP_0 = deltaFModified ; % -KModified_0*P0=0 (save calculations)
        fP = fP_0 ;
        err = inf ;
        NRIterNum = 0 ;
        while err > simParam.numSol.convergenceCriterion && ...
                NRIterNum <= simParam.numSol.maxIterationsNR
            
            cycInd = mod(round(cycIn/dCycGMS), lenCyc) + 1 ;
            invKModified = invKCycModified(:, :, cycInd) ; % negative inverse Jacobian
            dP = invKModified * fP ; % initial guess
            P = P + dP ;
            
            err = abs((fP' * dP) / (fP_0' * P)) ;
            
            cycIn = (u(6, ii) + P(6)) / (2*pi) ;
            cycInd = mod(round(cycIn/dCycGMS), lenCyc) + 1 ;
            KModified = KCycModified(:, :, cycInd) ;
            fP = deltaFModified - KModified * P + (KModified_0 - KModified) * u(:, ii) ;
            
            NRIterNum = NRIterNum + 1 ;
        end % of while
        
        %% Update the current step according to Newmark method %%
        du = P ;
        du_dot = (gamma/(beta*dt))*du - (gamma/beta)*u_dot(:,ii) + dt*(1-(gamma/(2*beta)))*u_ddot(:,ii) ;
        du_ddot = (1/(beta*(dt^2)))*du - (1/(beta*dt))*u_dot(:,ii) - (1/(2*beta))*u_ddot(:,ii) ;
        u(:, ii+1) = u(:, ii) + du ;
        u_dot(:, ii+1) = u_dot(:,ii) + du_dot ;
        u_ddot(:, ii+1) = u_ddot(:,ii) + du_ddot ;
        
        %% Avoid side effcets of significant errors %%
        cycIn = u(6, ii+1) / (2*pi);
        cycInd = mod(round(cycIn/dCycGMS), lenCyc) + 1 ;
        K = KCyc(:, :, cycInd) ;
        FexVctrT(:, ii+1) = M*u_ddot(:,ii+1) + C*u_dot(:,ii+1) + K*u(:,ii+1) ;
        
        %% Update Waitbar
        try if ~rem(ii, round(0.05*(lenT-1)))
                waitbar(ii / lenT, h, sprintf('Processing... %d%%', round(ii/lenT*100)));
            end % of if
        end % of try
    end % of for ii
    try close(h) ; end % close waitbar
    
    %% Trim the signals to discard the numerical oscillations untill convergence %%
    sigFirstInd = simParam.numSol.sigFirstInd ;
    u = u(:, sigFirstInd:end)' ;
    u_dot = u_dot(:, sigFirstInd:end)' ;
    u_ddot = u_ddot(:, sigFirstInd:end)' ;
    
end % of function 'numericalSolution'