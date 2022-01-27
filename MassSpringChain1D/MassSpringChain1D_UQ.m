%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Application of Chen and Horing method for material extrapolation
% From a potential to another potential
% Mass spring chain with N degrees of freedom
% |-----0-----0-----0----->
%       |->x1 |->x2 |->x3 |->lambda
% eta xdot = - DV + sqrt(2 eta kBT) Wdot
% force = DV
%
% Code 2 - Uncertainty quantification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

TClock_0 = clock;

%% Parameters setting for uncertainty quantification
PlotProp = 0;   % 0 for not to plot the linear variance propagation result 
                % 1 for plotting the linear variance propagation result 
                % Notice the linear variance propagation assumes starting
                % from Brownian motion (V1=0)

% V_1(u), vp_1: (potenial and pulling velocity) of simulated system
% V_2(u), vp_2: target system
% V_bias(u)= V_1 - V_2;

N=3;                       % Number of degress of freedom
eta=5;
kBT=1e-4;          
vp_1 = 0.01;                 % Pulling velocity for V1
vp_2 = 0.01;            	% Pulling velocity for V2

% Time discretization for Langevin simulation. 
NR = 1e3;             % Number of realizations
T = 1;               % Time interval [0,T]
dt = 1e-3;            % Time step for Langevin
Nt = T/dt;

% Definition of the potential
type = 'quartic';

if  strcmp(type,'quartic')
    % From quartic potential V1 to quartic potential V2
    % V(u) = 1/2 * k2 * u^2 + 1/4 * k4 * u^4
    k2_1=0; k4_1=00;          %   parameters for V1  
    k2_2=1; k4_2= 00;           %   parameters for V2
    DV1= @(x) DQuarticPotential(k2_1,k4_1,N,x);
    DV2= @(x) DQuarticPotential(k2_2,k4_2,N,x);
else
    disp('Wrong potential type');
    stop
end


% Load data file from Langevin simulation

DataFileName = ['_N' num2str(N,'%g') '_R' num2str(NR,'%g') '_k' num2str(k2_1,'%g') '_' num2str(k2_2,'%g') ...
    '_4k' num2str(k4_1,'%g') '_' num2str(k4_2,'%g') '_eta' num2str(eta,'%g') '_kT' num2str(kBT,'%g')...
    '_v' num2str(vp_1,'%g') '_' num2str(vp_2,'%g') ...
    '_dt' num2str(dt,'%g')];

load(['Result/DataLang' DataFileName '.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uncertainty quantification

% select the type of reference trajectory for potential expansion
type_xr = '0';                  % expand at 0 displacement
% type_xr = 'linear';           % expand at displacement field that is linear in space under given pulling volocity

Nt_est = 100;           % number of estimation time points, which is selected linearly within [0,T] at first
                        % and then updated later to finer the mesh at beginning time in log-scale
Nt_grid_est = 100;      % number of time discretization for UQ at given estimation time point t, 
                        % which uniformly discretize the time interval [0, t]
time_est = T/Nt_est:T/Nt_est:T; % mesh for estimation time points

% finer the mesh at beginning (in log-scale) if the current first point is too late
while time_est(1) >= 4 * dt
    Nt_est = Nt_est + 1;
    time_est = [time_est(1)/2, time_est];
end
if time_est(1) > 2 * dt
    Nt_est = Nt_est + 1;
    time_est = [2 * dt, time_est];
end

time_grid_est = zeros(Nt_grid_est+1,Nt_est);    % time_grid_est(:,iT) is the time grid that discretize [0, time_est(iT)]
for iT = 1:Nt_est
    time_grid_est(:,iT) = linspace(0, time_est(iT), Nt_grid_est+1)';
end
dt_est = time_est/Nt_grid_est;  % timestep for UQ at different estimation time point


% define the function for reference trajectory
if strcmp(type_xr, '0')
    fun_xr_t = @(t) zeros(N,length(t));
elseif strcmp(type_xr, 'linear')
    fun_xr_t = @(t) xr_t_linear(N, vp_1 * t);
else
    disp('Wrong type_xr!');
    stop
end

% define the functions for the gradient and double gradient of V_1, V_2
% and V_bias at reference trajectory x_r(t) (for nonlinear UQ) and 0 (for linear UQ)
if strcmp(type, 'quartic') 
    fun_DV1r_t = @(t) DV_x_Quartic(k2_1, k4_1, N, [fun_xr_t(t); vp_1*t]);
    fun_D2V1r_t = @(t) D2V_x_Quartic(k2_1, k4_1, N, [fun_xr_t(t); vp_1*t]);    
    
    fun_DV2r_t = @(t) DV_x_Quartic(k2_2, k4_2, N, [fun_xr_t(t); vp_2*t]);
    fun_D2V2r_t = @(t) D2V_x_Quartic(k2_2, k4_2, N, [fun_xr_t(t); vp_2*t]);
    
    fun_DVbr_t = @(t) fun_DV1r_t(t) - fun_DV2r_t(t);
    fun_D2Vbr_t = @(t) fun_D2V1r_t(t) - fun_D2V2r_t(t);
    
    
    fun_DV10_t = @(t) DV_x_Quartic(k2_1, k4_1, N, [0*fun_xr_t(t); vp_1*t]);
    fun_D2V10_t = @(t) D2V_x_Quartic(k2_1, k4_1, N, [0*fun_xr_t(t); vp_1*t]);    
    
    fun_DV20_t = @(t) DV_x_Quartic(k2_2, k4_2, N, [0*fun_xr_t(t); vp_2*t]);
    fun_D2V20_t = @(t) D2V_x_Quartic(k2_2, k4_2, N, [0*fun_xr_t(t); vp_2*t]);
    
    fun_DVb0_t = @(t) fun_DV10_t(t) - fun_DV20_t(t);
    fun_D2Vb0_t = @(t) fun_D2V10_t(t) - fun_D2V20_t(t);
    
else 
    disp('Wrong potential type!');
    stop
end

% Calculate A, b, c factors for mean (with label y_t) and variance (with label vy_t) 
% of P_bias at all estimation time points
% For each estimation time point, A, b, c are matrix, vector and scalar, respectively 
detAy_t = ones(1,Nt_est);       % determinant of Ay_t
detAvy_t = ones(1,Nt_est);      % determinent of Avy_t
cy_t = zeros(1,Nt_est);
cvy_t = zeros(1,Nt_est);
b_A_inv_b_y_t = zeros(1,Nt_est);    % b^T A^{-1} b for label y_t 
bv_Av_inv_bv_y_t = zeros(1,Nt_est); % b^T A^{-1} b for label vy_t 

for iT = 1:Nt_est    
    xrt = fun_xr_t(time_grid_est(:,iT)');   % an N * (Nt_grid_est+1) matrix for x_r(t)
    dxrt = diff(xrt, 1, 2) / dt_est(iT);    % an N * Nt_grid_est matrix for dx_r(t)/dt
    
    dV1rt = fun_DV1r_t(time_grid_est(:,iT)');   % an N * (Nt_grid_est+1) matrix for \nabla V1
    d2V1rt = fun_D2V1r_t(time_grid_est(:,iT)'); % an N*N * (Nt_grid_est+1) matrix for \nabla \nabla V1
    dV2rt = fun_DV2r_t(time_grid_est(:,iT)');   % an N * (Nt_grid_est+1) matrix for \nabla V2
    d2V2rt = fun_D2V2r_t(time_grid_est(:,iT)'); % an N*N * (Nt_grid_est+1) matrix for \nabla \nabla V2
    
    % initialize and calculate A, b, c
    by = zeros(N*Nt_grid_est, 1);
    bvy = zeros(N*Nt_grid_est, 1);
    Ay = zeros(N*Nt_grid_est, N*Nt_grid_est);
    Avy = zeros(N*Nt_grid_est, N*Nt_grid_est);
    
    cy_t(iT) = - 1/(2* 2*kBT*eta) * sum(sum( ( eta * dxrt + dV2rt(:,1:end-1) ).^2)) * dt_est(iT);
    cvy_t(iT) = 2*cy_t(iT) + 1/(2* 2*kBT*eta) * sum(sum( ( eta * dxrt + dV1rt(:,1:end-1) ).^2)) * dt_est(iT);
    
    for nT = 1:Nt_grid_est-1
        by(((nT-1)*N+1):(nT*N)) = dt_est(iT) * ( - (dV2rt(:,nT+1) - dV2rt(:,nT)) / dt_est(iT) ...
            - eta * (dxrt(:,nT+1) - dxrt(:,nT)) / dt_est(iT) ...
            + d2V2rt(:,:,nT+1) * dxrt(:,nT+1) + 1/eta * d2V2rt(:,:,nT+1) * dV2rt(:,nT+1) );
        
        bvy(((nT-1)*N+1):(nT*N)) = 2* by(((nT-1)*N+1):(nT*N)) ...
            - dt_est(iT) * ( - (dV1rt(:,nT+1) - dV1rt(:,nT)) / dt_est(iT) ...
            - eta * (dxrt(:,nT+1) - dxrt(:,nT)) / dt_est(iT) ...
            + d2V1rt(:,:,nT+1) * dxrt(:,nT+1) + 1/eta * d2V1rt(:,:,nT+1) * dV1rt(:,nT+1) );
    end
    by(((Nt_grid_est-1)*N+1):(Nt_grid_est*N)) = eta * dxrt(:,Nt_grid_est) + dV2rt(:,Nt_grid_est);
    bvy(((Nt_grid_est-1)*N+1):(Nt_grid_est*N)) = 2 * by(((Nt_grid_est-1)*N+1):(Nt_grid_est*N)) ...
        - ( eta * dxrt(:,Nt_grid_est) + dV1rt(:,Nt_grid_est) );
    
    by = -sqrt(dt_est(iT)/(2*kBT*eta)) * by;
    bvy = -sqrt(dt_est(iT)/(2*kBT*eta)) * bvy;

    
    for nT = 1:Nt_grid_est-1
        Gamma1 = eye(N) - dt_est(iT)/eta * d2V1rt(:,:,nT+1);
        Gamma2 = eye(N) - dt_est(iT)/eta * d2V2rt(:,:,nT+1);
        Gamma1sq = Gamma1 * Gamma1;
        Gamma2sq = Gamma2 * Gamma2;
        
        Ay( ((nT-1)*N+1):(nT*N), ((nT-1)*N+1):(nT*N) ) = eye(N) + Gamma2sq;
        Ay( ((nT-1)*N+1):(nT*N), ((nT)*N+1):((nT+1)*N) ) = - Gamma2;
        Ay( ((nT)*N+1):((nT+1)*N), ((nT-1)*N+1):(nT*N) ) = - Gamma2;     
        
        Avy( ((nT-1)*N+1):(nT*N), ((nT-1)*N+1):(nT*N) ) = 2 * (eye(N) + Gamma2sq) - (eye(N) + Gamma1sq);
        Avy( ((nT-1)*N+1):(nT*N), ((nT)*N+1):((nT+1)*N) ) = - 2 * Gamma2 + Gamma1;
        Avy( ((nT)*N+1):((nT+1)*N), ((nT-1)*N+1):(nT*N) ) = - 2 * Gamma2 + Gamma1;  
    end
    Ay( ((Nt_grid_est-1)*N+1):(Nt_grid_est*N), ((Nt_grid_est-1)*N+1):(Nt_grid_est*N) ) = eye(N);
    Avy( ((Nt_grid_est-1)*N+1):(Nt_grid_est*N), ((Nt_grid_est-1)*N+1):(Nt_grid_est*N) ) = eye(N);
        
    detAy_t(iT) = det(Ay);
    detAvy_t(iT) = det(Avy);
    
    b_A_inv_b_y_t(iT) = by' * (Ay\by);
    bv_Av_inv_bv_y_t(iT) = bvy' * (Avy\bvy);
    
end

% Estimated mean and variance of P_bias based on nonlinear UQ (_2nd indicates quadratic expansion)
% Notice that mean should be 1
mean_Pfactor_2nd = 1./sqrt(detAy_t) .* exp(b_A_inv_b_y_t/2 + cy_t);
var_Pfactor_2nd = 1./sqrt(detAvy_t) .* exp(bv_Av_inv_bv_y_t/2 + cvy_t) - mean_Pfactor_2nd.^2;


if PlotProp == 1
    % Calculate the variance of P_bias based on linear UQ
    % Notice that the time discretization for linear UQ is different from nonlinear UQ
    % Linear UQ discretize [0,T] into Nt_pre sections and use dt_pre as timestep for 
    % all estimation time points dt_pre*(1:Nt_pre)
    
    dt_pre = 2*dt;                  % timestep for linear UQ (same for all estimation time points)
    Nt_pre = Nt * (dt/dt_pre);
    time_pre = dt_pre:dt_pre:T;     % estimation time points for linear UQ

    if strcmp(type, 'quartic')
        if k2_1~=0 || k4_1~=0       % this code only works for V_1 = 0
            disp('V_1 should be 0 for linear UQ!');
            stop
        end
    end

    dVb0t_pre = fun_DVb0_t(time_pre-dt_pre);
    d2Vb0t_pre = fun_D2Vb0_t(time_pre-dt_pre);

    g_num = 0*time_pre;         % g_num = sum_ti dVb0(ti).^2
    g_num(1) = dVb0t_pre(1).^2 * dt_pre;
    for iT = 2:Nt_pre
        g_num(iT) = g_num(iT-1) + dVb0t_pre(iT).^2  * dt_pre;
    end

    h_num = zeros(N, Nt_pre,Nt_pre);    % h_num(:,i,j) = sum_{t=i}^{t=j} dVb0(t) \dot d2Vb0(t) dt_pre
                                        % Notice that h_num(:,i,j) = 0 if i > j
    dVb0t_d2V0t_pre = zeros(N,Nt_pre);  % inner product dVb0(t) \dot d2Vb0(t) for all t
    for jT = 1:Nt_pre
        dVb0t_d2V0t_pre(:,jT) = dVb0t_pre(:,jT)' * d2Vb0t_pre(:,:,jT);
    end
    for iT = 1:Nt_pre
        for jT=iT+1:Nt_pre
            h_num(:,iT,jT) = h_num(:,iT,jT-1) + dVb0t_d2V0t_pre(:,jT) * dt_pre;
        end
    end

    g2_num = zeros(1,Nt_pre);
    for iT = 1:Nt_pre
        g2_num(iT) = sum(sum((dVb0t_pre(:,1:iT) - 1/eta * h_num(:,1:iT,iT)).^2)) * dt_pre;
        % two sum() indicate the summations over time interval [0, iT*dt_pre] and N particles
    end
    
    % estimated variance of P_bias using linear UQ (_prop for linear variance propagation)
    var_Pfactor_prop = 1/(2*kBT*eta) * exp(-(2/4*kBT*eta) * g_num) .* g2_num;
    % where the exponential term is (P_bias)^2 when all random noises are 0
end


% Check whether the nonlinear UQ mean of P_bias is 1
figure
plot(time_est,mean_Pfactor_2nd,'LineWidth',3,'DisplayName','2nd')
hold off
set(gca,'yscale','log','xscale','log')
legend('Location','NorthWest')
legend('Interpreter','latex')
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$<\mathcal{P}_b>$','Interpreter','latex','FontSize',20)


% Plot the variance of P_bias from Langevin data, nonlinear UQ and linear UQ (if PlotProp==1)
figure
plot(time,var_Pfactor,'LineWidth',3,'DisplayName','data')
hold on
plot(time_est, var_Pfactor_2nd, '--','LineWidth',3,'DisplayName','nonlinear')
hold on
if PlotProp == 1
    plot(time_pre, var_Pfactor_prop, '--','LineWidth',3,'DisplayName','linear')
end
hold off
set(gca,'FontSize',20,'xscale','log','yscale','log');
legend('Location','NorthWest')
legend('Interpreter','latex','FontSize',15)
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$Var(e^{-\beta I_b})$','Interpreter','latex','FontSize',20)




Nfactor = normalization/NR;         % the finite average of P_bias
err_Nfactor = abs(Nfactor - 1);     % the deviation of Nfactor from 1
err_Nfactor_2nd = sqrt(var_Pfactor_2nd/NR);         % the estimated deviation (standar error) of Nfactor based on nonlinear UQ
if PlotProp == 1
    err_Nfactor_prop = sqrt(var_Pfactor_prop/NR);   % the estimated deviation (standar error) of Nfactor based on linear UQ
end


figure
plot(time,Nfactor,'LineWidth',3)
set(gca,'FontSize',20, 'YScale', 'log');
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$\mathcal{N}$','Interpreter','latex','FontSize',20)


figure
plot(time,err_Nfactor,'LineWidth',3,'DisplayName','data')
hold on
plot(time_est,err_Nfactor_2nd, '--', 'LineWidth',3,'DisplayName','nonlinear')
hold on
if PlotProp == 1
    plot(time_pre,err_Nfactor_prop, '--', 'LineWidth',3,'DisplayName','linear')
end
hold off
set(gca,'FontSize',15, 'YScale', 'log');
legend('Location','SouthEast','NumColumns',2)
legend('Interpreter','latex')
axis([-inf, inf, -inf, 2])
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$|\mathcal{N}-1|$','Interpreter','latex','FontSize',20)



% Save data

EstFileName = ['_N' num2str(N,'%g') '_R' num2str(NR,'%g') '_k' num2str(k2_1,'%g') '_' num2str(k2_2,'%g') ...
    '_4k' num2str(k4_1,'%g') '_' num2str(k4_2,'%g') '_eta' num2str(eta,'%g') '_kT' num2str(kBT,'%g')...
    '_v' num2str(vp_1,'%g') '_' num2str(vp_2,'%g') ...
    '_dt' num2str(dt,'%g') '_NG' num2str(Nt_grid_est, '%g') '_xr' type_xr(1)];

save(['Result/DataEst' EstFileName '.mat'], 'dt_est', 'Nt_est', 'Nt_grid_est', 'time_est', 'time_grid_est',  ...
    'cy_t', 'cvy_t', 'detAvy_t', 'detAy_t', 'b_A_inv_b_y_t', 'bv_Av_inv_bv_y_t', ...
    'type_xr', 'xrt', 'mean_Pfactor_2nd', 'var_Pfactor_2nd',...
    'Nfactor', 'err_Nfactor', 'err_Nfactor_2nd')
if PlotProp == 1
    save(['Result/DataEst' EstFileName '.mat'], 'dt_pre', 'Nt_pre', 'time_pre', 'var_Pfactor_prop', 'err_Nfactor_prop', '-append')
end


FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
mkdir('Result')
mkdir('Result/figures')


savefig(FigList,['Result/figures/FigEst' EstFileName '.fig'])


TClock_f = clock;
calTime_Est = etime(TClock_f,TClock_0)

save(['Result/DataEst' EstFileName '.mat'], 'calTime_Est', '-append')