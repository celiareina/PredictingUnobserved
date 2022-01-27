%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Application of Chen and Horing method for material extrapolation
% From a potential V1 to another potential V2
% 2D box with size L, periodic BC and N particles with A, B two specieses.
% eta rdot = - DV + sqrt(2 eta kBT) Wdot
% force = DV
%
% Code 2 - Uncertainty quantification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

TClock_0 = clock;


% type of potential
type = 'Hertzian';
% type = 'LJ';


% Parameters of the system
% _1: system we simulate (only one)
% _2: systems we want to predict (multiple, number NV)

if  strcmp(type,'Hertzian')
    % From Hertzian to Hertzian    
    % V(r) = 0.4 * Vm * (max((1 - |r|/sigma), 0))^2.5
    
    FolderName = 'Result/2D_Hertzian/';     % Data folder name
    
    % Data folder & file name & potential amplitude for V1
    InputFileName = 'hertz.mat';
    Vm_1 = 10;

%     InputFileName = 'hertz_lowt.mat';
%     Vm_1 = 1;

    load([FolderName InputFileName]);   % load data (including system parameters & variance of P_bias)
    
    rad_A = 0.5;    %  radius of particle A
    rad_B = 0.7;    % radius of particle B
    index_A = 1:5;	% indices for particle A

elseif strcmp(type, 'LJ')
    % From Lennard-Jones (LJ) to LJ    
    % V(r) = 4 * Vm * ( (sigma / |r|)^12 - (sigma / |r|)^6 )
    
    FolderName = 'Result/2D_LJ/';       % Data folder name
    
	% Data folder & file name & potential amplitude between A-A for V1 
    InputFileName = 'lj.mat';
    Vm_AA_1 = 0.1;

%     InputFileName = 'lj_lowt.mat';
%     Vm_AA_1 = 0.01;

    load([FolderName InputFileName]);
    
    index_A = 5:10;     % index for particle A
    
    % interatomic length (sigma) and amplitudes (Vm) between different species (A/B)
    sigma_AA = 1;	
    sigma_AB = 0.8 * sigma_AA;
    sigma_BB = 0.88 * sigma_AA;
    Vm_AB_1 = 1.5 * Vm_AA_1;
    Vm_BB_1 = 0.5 * Vm_AA_1;
    
    sigma = [sigma_AA; sigma_AB; sigma_BB];
    Vm_1 = [Vm_AA_1; Vm_AB_1; Vm_BB_1];
end

% Convert the parameter names from the loading file
dVms_all = dVs;         % potential amplitude difference Vm_2 - Vm_1, size of NV
Pb_mean_all = Ib_mean;  % mean of P_bias for all different dVms, size of [NV, Nt]
Pb_var_all = Ib_var;    % variance of P_bias for all different dVms, size of [NV, Nt]
NDim = double(dim);     % dimension of problem, which is not restricted as 2
NR = double(N);         % number of realization
N = double(num_parts);  % number of particles
eta = visc;             % viscosity
T = time(end);          % total time
kBT = 1/beta;           % temperature


[NV, Nt] = size(Pb_mean_all);   % NV as the number of dVms (or target systems), Nt as the number of timesteps
NDoF = N * NDim;        % number of degree of freedom (DoF)

init_pos = init_pos(:,1:NDim);  % initial position for all particles in all directions, size of [N, NDim] with order
                                % [x1, y1, ...; x2, y2, ...; ...; xN, yN, ...]
x_ini = init_pos';
x_ini = x_ini(:);               % convert to 1D array with size NDoF and order
                                % [x1, y1, ..., x2, y2, ..., ..., xN, yN, ...]

Vm_2_all = Vm_1 * (1 - dVms_all);   % all the target potential amplitudes

clear Ib_mean Ib_var dVs dim num_parts visc


iV_index = [2, NV];         % select the index/indices of target systems/potentials for applying UQ 
NumV = length(iV_index);    % number of selected target potentials
% do loop for all selected target system
for iV = 1:NumV
    close all;      % close the plotted figures for previous target system
    TClock_1 = clock;
    
    % the target potential and corresponding mean and variance of P_bias
    Vm_2 = Vm_2_all(:,iV_index(iV));
    Pb_mean = Pb_mean_all(iV_index(iV),:);
    Pb_var = Pb_var_all(iV_index(iV),:);

    
    % Plot the initial position with different colors for specieses A & B
    figure
    for i = 1:N
        if ismember(i, index_A)
            scatter(init_pos(i,1), init_pos(i,2), 'filled', 'DisplayName',num2str(i))
        else
            scatter(init_pos(i,1), init_pos(i,2), 'DisplayName',num2str(i))
        end
        hold on
    end
    hold off
    set(gca,'FontSize',15);
    legend('Location','SouthEast','NumColumns',2)
    legend('Interpreter','latex')
    axis([-L/2, L/2, -L/2, L/2])
    xlabel('$x$','Interpreter','latex','FontSize',20)
    ylabel('$y$','Interpreter','latex','FontSize',20)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% estimation of normalization factor

    % select the type of reference trajectory for potential expansion
    type_xr = '0';          % expand at 0 displacement
    
    
    Nt_est = 100;           % number of estimation time points, which is selected linearly within [0,T] at first
                            % and then updated later to finer the mesh at beginning time in log-scale
    Nt_grid_est = 100;     	% number of time discretization for UQ at given estimation time point t, 
                            % which uniformly discretize the time interval [0, t]
    time_est = (T/Nt_est:T/Nt_est:T);   % mesh for estimation time points
    
    % finer the mesh at beginning (in log-scale) if the current first point is too late
    while time_est(1) >= 2*time(1)
        Nt_est = Nt_est + 1;
        time_est = [time_est(1)/2, time_est];
    end
    if time_est(1) > time(1)
        Nt_est = Nt_est + 1;
        time_est = [time(1), time_est];
    end
    
    time_grid_est = zeros(Nt_grid_est+1,Nt_est);    % time_grid_est(:,iT) is the time grid that discretizes [0, time_est(iT)]
    for iT = 1:Nt_est
        time_grid_est(:,iT) = linspace(0, time_est(iT), Nt_grid_est+1)';
    end
    dt_est = time_est/Nt_grid_est;      % timestep for UQ at different estimation time point


    % define the function for reference trajectory
    if strcmp(type_xr, '0')
        fun_xr_t = @(t) zeros(NDoF,length(t)) + x_ini;
    else
        disp('Wrong type_xr!');
        stop
    end

    % define the functions for the gradient and double gradient of V_1, V_2
    % and V_bias at reference trajectory x_r(t)
    if strcmp(type, 'Hertzian')
        fun_DV1r_t = @(t) DV_x_Hertzian(Vm_1, rad_A, rad_B, N, index_A, NDim, fun_xr_t(t), L);
        fun_D2V1r_t = @(t) D2V_x_Hertzian(Vm_1, rad_A, rad_B, N, index_A, NDim, fun_xr_t(t), L);    

        fun_DV2r_t = @(t) DV_x_Hertzian(Vm_2, rad_A, rad_B, N, index_A, NDim, fun_xr_t(t), L);
        fun_D2V2r_t = @(t) D2V_x_Hertzian(Vm_2, rad_A, rad_B, N, index_A, NDim, fun_xr_t(t), L);

        fun_DVbr_t = @(t) fun_DV1r_t(t) - fun_DV2r_t(t);
        fun_D2Vbr_t = @(t) fun_D2V1r_t(t) - fun_D2V2r_t(t);
        
    elseif strcmp(type, 'LJ')
        fun_DV1r_t = @(t) DV_x_LJ(Vm_1, sigma, N, index_A, NDim, fun_xr_t(t), L);
        fun_D2V1r_t = @(t) D2V_x_LJ(Vm_1, sigma, N, index_A, NDim, fun_xr_t(t), L);    

        fun_DV2r_t = @(t) DV_x_LJ(Vm_2, sigma, N, index_A, NDim, fun_xr_t(t), L);
        fun_D2V2r_t = @(t) D2V_x_LJ(Vm_2, sigma, N, index_A, NDim, fun_xr_t(t), L);

        fun_DVbr_t = @(t) fun_DV1r_t(t) - fun_DV2r_t(t);
        fun_D2Vbr_t = @(t) fun_D2V1r_t(t) - fun_D2V2r_t(t);
        
    else
        disp('Wrong potential type!');
        stop
    end

    % Calculate A, b, c factors for mean (with label y_t) and variance (with label vy_t) 
    % of P_bias at all estimation time points
    % For each estimation time point, A, b, c are matrix, vector and scalar, respectively 
    detAy_t = ones(1,Nt_est);   % determinant of Ay_t
    detAvy_t = ones(1,Nt_est);  % determinent of Avy_t
    cy_t = zeros(1,Nt_est);
    cvy_t = zeros(1,Nt_est);
    b_A_inv_b_y_t = zeros(1,Nt_est);    % b^T A^{-1} b for label y_t 
    bv_Av_inv_bv_y_t = zeros(1,Nt_est); % b^T A^{-1} b for label vy_t 
    
    for iT = 1:Nt_est
        xrt = fun_xr_t(time_grid_est(:,iT)');	% an N * (Nt_grid_est+1) matrix for x_r(t)
        dxrt = diff(xrt, 1, 2) / dt_est(iT);    % an N * Nt_grid_est matrix for dx_r(t)/dt

        dV1rt = fun_DV1r_t(time_grid_est(:,iT)');   % an NDoF * (Nt_grid_est+1) matrix for \nabla V1
        d2V1rt = fun_D2V1r_t(time_grid_est(:,iT)'); % an NDoF*NDoF * (Nt_grid_est+1) matrix for \nabla \nabla V1
        dV2rt = fun_DV2r_t(time_grid_est(:,iT)');   % an NDoF * (Nt_grid_est+1) matrix for \nabla V2
        d2V2rt = fun_D2V2r_t(time_grid_est(:,iT)'); % an NDoF*NDoF * (Nt_grid_est+1) matrix for \nabla \nabla V2

        % initialize and calculate A, b, c
        by = zeros(NDoF*Nt_grid_est, 1);
        bvy = zeros(NDoF*Nt_grid_est, 1);
        Ay = zeros(NDoF*Nt_grid_est, NDoF*Nt_grid_est);
        Avy = zeros(NDoF*Nt_grid_est, NDoF*Nt_grid_est);

        cy_t(iT) = - 1/(2* 2*kBT*eta) * sum(sum( ( eta * dxrt + dV2rt(:,1:end-1) ).^2)) * dt_est(iT);
        cvy_t(iT) = 2*cy_t(iT) + 1/(2* 2*kBT*eta) * sum(sum( ( eta * dxrt + dV1rt(:,1:end-1) ).^2)) * dt_est(iT);

        for nT = 1:Nt_grid_est-1
            by(((nT-1)*NDoF+1):(nT*NDoF)) = dt_est(iT) * ( - (dV2rt(:,nT+1) - dV2rt(:,nT)) / dt_est(iT) ...
                - eta * (dxrt(:,nT+1) - dxrt(:,nT)) / dt_est(iT) ...
                + d2V2rt(:,:,nT+1) * dxrt(:,nT+1) + 1/eta * d2V2rt(:,:,nT+1) * dV2rt(:,nT+1) );

            bvy(((nT-1)*NDoF+1):(nT*NDoF)) = 2* by(((nT-1)*NDoF+1):(nT*NDoF)) ...
                - dt_est(iT) * ( - (dV1rt(:,nT+1) - dV1rt(:,nT)) / dt_est(iT) ...
                - eta * (dxrt(:,nT+1) - dxrt(:,nT)) / dt_est(iT) ...
                + d2V1rt(:,:,nT+1) * dxrt(:,nT+1) + 1/eta * d2V1rt(:,:,nT+1) * dV1rt(:,nT+1) );
        end
        by(((Nt_grid_est-1)*NDoF+1):(Nt_grid_est*NDoF)) = eta * dxrt(:,Nt_grid_est) + dV2rt(:,Nt_grid_est);
        bvy(((Nt_grid_est-1)*NDoF+1):(Nt_grid_est*NDoF)) = 2 * by(((Nt_grid_est-1)*NDoF+1):(Nt_grid_est*NDoF)) ...
            - ( eta * dxrt(:,Nt_grid_est) + dV1rt(:,Nt_grid_est) );

        by = -sqrt(dt_est(iT)/(2*kBT*eta)) * by;
        bvy = -sqrt(dt_est(iT)/(2*kBT*eta)) * bvy;


        for nT = 1:Nt_grid_est-1
            Gamma1 = eye(NDoF) - dt_est(iT)/eta * d2V1rt(:,:,nT+1);
            Gamma2 = eye(NDoF) - dt_est(iT)/eta * d2V2rt(:,:,nT+1);
            Gamma1sq = Gamma1 * Gamma1;
            Gamma2sq = Gamma2 * Gamma2;

            Ay( ((nT-1)*NDoF+1):(nT*NDoF), ((nT-1)*NDoF+1):(nT*NDoF) ) = eye(NDoF) + Gamma2sq;
            Ay( ((nT-1)*NDoF+1):(nT*NDoF), ((nT)*NDoF+1):((nT+1)*NDoF) ) = - Gamma2;
            Ay( ((nT)*NDoF+1):((nT+1)*NDoF), ((nT-1)*NDoF+1):(nT*NDoF) ) = - Gamma2;     

            Avy( ((nT-1)*NDoF+1):(nT*NDoF), ((nT-1)*NDoF+1):(nT*NDoF) ) = 2 * (eye(NDoF) + Gamma2sq) - (eye(NDoF) + Gamma1sq);
            Avy( ((nT-1)*NDoF+1):(nT*NDoF), ((nT)*NDoF+1):((nT+1)*NDoF) ) = - 2 * Gamma2 + Gamma1;
            Avy( ((nT)*NDoF+1):((nT+1)*NDoF), ((nT-1)*NDoF+1):(nT*NDoF) ) = - 2 * Gamma2 + Gamma1;  
        end
        Ay( ((Nt_grid_est-1)*NDoF+1):(Nt_grid_est*NDoF), ((Nt_grid_est-1)*NDoF+1):(Nt_grid_est*NDoF) ) = eye(NDoF);
        Avy( ((Nt_grid_est-1)*NDoF+1):(Nt_grid_est*NDoF), ((Nt_grid_est-1)*NDoF+1):(Nt_grid_est*NDoF) ) = eye(NDoF);

        detAy_t(iT) = det(Ay);
        detAvy_t(iT) = det(Avy);

        b_A_inv_b_y_t(iT) = by' * (Ay\by);
        bv_Av_inv_bv_y_t(iT) = bvy' * (Avy\bvy);
        
    end

    % Estimated mean and variance of P_bias based on nonlinear UQ (_2nd indicates quadratic expansion)
    % Notice that mean should be 1
    mean_Pfactor_2nd = exp(b_A_inv_b_y_t/2 + cy_t);
    var_Pfactor_2nd = 1./sqrt(detAvy_t) .* exp(bv_Av_inv_bv_y_t/2 + cvy_t) - mean_Pfactor_2nd.^2;

    
    % Check whether the nonlinear UQ mean of P_bias is 1
    figure
    plot(time_est,mean_Pfactor_2nd,'LineWidth',3,'DisplayName','2nd')
    hold off
    set(gca,'yscale','log','xscale','log')
    legend('Location','NorthWest')
    legend('Interpreter','latex')
    xlabel('$t$','Interpreter','latex','FontSize',20)
    ylabel('$<\mathcal{P}_b>$','Interpreter','latex','FontSize',20)

    
    % Plot the variance of P_bias from Langevin data and nonlinear UQ
    figure
    plot(time,Pb_var,'LineWidth',3,'DisplayName','data')
    hold on
    plot(time_est, var_Pfactor_2nd, '--','LineWidth',3,'DisplayName','prediction')
    hold off
    set(gca,'FontSize',20,'xscale','log','yscale','log');
    legend('Location','NorthWest')
    legend('Interpreter','latex','FontSize',15)
    xlabel('$t$','Interpreter','latex','FontSize',20)
    ylabel('$Var(e^{-\beta I_b})$','Interpreter','latex','FontSize',20)
    
    
    std_Pfactor_2nd = sqrt(var_Pfactor_2nd);
    for iT = 1:Nt_est
        if abs(real(var_Pfactor_2nd(iT))) < abs(5*imag(var_Pfactor_2nd(iT)))
            % when the predicted standard deviation from UQ is not real number (variance is negative and not relialbe),
            % set the std from this estimation time points to the end as infinity for not showing them in the plot
            std_Pfactor_2nd(iT:end) = Inf;
            break;
        end
    end
    
    figure
    plot(time,sqrt(Pb_var),'LineWidth',3,'DisplayName','data')
    hold on
    plot(time_est, std_Pfactor_2nd, '--','LineWidth',3,'DisplayName','prediction')
    hold off
    set(gca,'FontSize',20,'xscale','log','yscale','log');
    legend('Location','NorthWest')
    legend('Interpreter','latex','FontSize',15)
    xlabel('$t$','Interpreter','latex','FontSize',20)
    ylabel('$\sigma_{\mathcal{P}_{bias}}$','Interpreter','latex','FontSize',20)


    Nfactor = Pb_mean;                      % the observed finite average of P_bias
    err_Nfactor = abs(Nfactor - 1);         % the observed deviation of Nfactor from 1
%     err_Nfactor_2nd = sqrt(var_Pfactor_2nd/double(NR));
    err_Nfactor_2nd = std_Pfactor_2nd/sqrt(double(NR));     % the estimated deviation (standar error) of Nfactor based on UQ

    
    figure
    plot(time,Nfactor,'LineWidth',3)
    set(gca,'FontSize',20, 'YScale', 'log');
    xlabel('$t$','Interpreter','latex','FontSize',20)
    ylabel('$\mathcal{N}$','Interpreter','latex','FontSize',20)


    figure
    plot(time,err_Nfactor,'LineWidth',3,'DisplayName','data')
    hold on
    plot(time_est,err_Nfactor_2nd, '--', 'LineWidth',3,'DisplayName','prediction')
    hold off
    set(gca,'FontSize',15, 'YScale', 'log');
    legend('Location','SouthEast','NumColumns',2)
    legend('Interpreter','latex')
    xlabel('$t$','Interpreter','latex','FontSize',20)
    ylabel('$|\mathcal{N}-1|$','Interpreter','latex','FontSize',20)


    


    if strcmp(type, 'Hertzian')
        EstFileName = ['_N' num2str(N,'%g') '_' num2str(NDim) 'DH' '_R' num2str(NR,'%g')...
            '_Vm' num2str(Vm_1,'%g') '_' num2str(Vm_2, '%g')...
            '_eta' num2str(eta,'%g') '_kT' num2str(kBT,'%g')...
            '_T' num2str(T,'%g') '_NG' num2str(Nt_grid_est, '%g') '_xr' type_xr(1)];
    elseif strcmp(type, 'LJ')
        EstFileName = ['_N' num2str(N,'%g') '_' num2str(NDim) 'DLJ' '_R' num2str(NR,'%g')...
            '_Vm' num2str(Vm_1(1),'%g') '_' num2str(Vm_2(1), '%g')...
            '_eta' num2str(eta,'%g') '_kT' num2str(kBT,'%g')...
            '_T' num2str(T,'%g') '_NG' num2str(Nt_grid_est, '%g') '_xr' type_xr(1)];
    else
        disp('Wrong potential type!');
        stop
    end

    save([FolderName 'DataEst' EstFileName '.mat'], 'dt_est', 'Nt_est', 'Nt_grid_est', 'time_est', 'time_grid_est',  ...
        'cy_t', 'cvy_t', 'detAvy_t', 'detAy_t', 'b_A_inv_b_y_t', 'bv_Av_inv_bv_y_t', ...
        'type_xr', 'xrt', 'mean_Pfactor_2nd', 'var_Pfactor_2nd',...
        'Nfactor', 'err_Nfactor', 'err_Nfactor_2nd')



    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

    savefig(FigList,[FolderName 'FigEst' EstFileName '.fig'])   % save all the plotted figures for current target system


    TClock_f = clock;
    calTime_Est = etime(TClock_f,TClock_1)
    save([FolderName 'DataEst' EstFileName '.mat'], 'calTime_Est', '-append')


end

