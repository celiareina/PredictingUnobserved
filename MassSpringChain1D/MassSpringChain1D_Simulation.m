%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predict non-equilibirum behavior for target system/process 2 from
% simulated system/process 1.
% 1D mass spring chain with N degrees of freedom
% |-----0-----0-----0----->
%       |->x1 |->x2 |->x3 |->lambda
% eta \dot{x} = - DV + sqrt(2 eta kBT) \dot{\xi}
% force = DV
% lambda = vp * t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

TClock_0 = clock;

%% Direct Langevin simulations of a mass spring chain with N degrees of freedom and quadratic potentials

%Parameters of the system
% V_1(u), vp_1: (potenial and pulling velocity) of simulated system
% V_2(u), vp_2: target system
% V_bias(u)= V_1 - V_2;
N=3;                       % Number of degress of freedom
eta=5;                      % viscosity
kBT=1e-4;                   % temperature
vp_1 = 0.01;                % Pulling velocity for simulated system
vp_2 = 0.01;                % Pulling velocity for target system

% Time discretization for Langevin simulation
NR = 1e3;             % Number of realizations
T = 1;               % Total time for process
dt = 1e-3;            % Time step for Langevin
Nt = T/dt;

% Definition of the potential type
type = 'quartic';

% Whether to save factor P_bias for all realization in output file
PfactorOutput = 1;    % 0 for not saving P_bias (Pfactor)
                      % 1 for saving P_bias in file

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Langevin simulation for V1, V2 and prediction for V2 from V1
% _2 if V2
% _1 if V1
% _bias for predicting V2 from V1
% Wave work average
% Force external average force
% r position of all particles with time

r_1=zeros(N,Nt);                % Position of all particle
r_2=zeros(N,Nt);
r_bias=zeros(N,Nt);
Wave_2 = zeros(1,Nt);           % Work average
Wave_1 = zeros(1,Nt);
Wave_bias = zeros(1,Nt);
Force_2 = zeros(1,Nt);          % External force
Force_1 = zeros(1,Nt);          
Force_bias = zeros(1,Nt); 

factor=zeros(NR,Nt+1);          % The integrand factor I_bias
Pfactor=zeros(NR,Nt+1);         % Factor P_bias = exp(-beta I_bias)
normalization=zeros(1,Nt);      % Normalization factor, the summation of P_bias over all realizations


% Start Langevin simulation
r0 = zeros(N,1);        % Initial condition
for r=1:NR     % going over realizations
    
   % Integrating in time with target V_2 and vp_2
   % for validation purposes
   
   % Initialization
   rold=[r0;0];         % adding lambda(0)=0 at the end
   W = zeros(1,Nt+1);   % Work
   
   for j=1:Nt   % going over time
       lambda_2=vp_2*dt*j;      % displacement for external protocal 
       rold(N+1,1)=lambda_2;    % setting lambda(j)=lambda_2 at the end
       dnoise=sqrt(2*eta*kBT*dt)*randn(N,1);    % Gaussian noise
       DV=DV2(rold);            
       force=DV(1:N);       % potential gradient / force on all particles, positive when pointing leftward
       force_end=DV(N+1);   % external force, positive when pointing rightward
       
       % Euler-Maruyama
       rnew(1:N,1) = rold(1:N) - force*dt/eta + dnoise/eta;
       
       % Calculating observables
       W(j+1)= W(j) + force_end*vp_2*dt;
       Force_2(j) = Force_2(j) + force_end;
       Wave_2(j) = Wave_2(j) + W(j+1);
       r_2(:,j)=r_2(:,j) + rnew(1:N);
       
       %update
       rold=[rnew;0];
   end
   
   
   % Integrating in time with simulated V_1 and vp_1 accounting for bias
   
   % Initialization
   rold_1=[r0;0];   
   rold_2=[r0;0];   
   W_1 = zeros(1,Nt+1);
   W_2 = zeros(1,Nt+1);

   for j=1:Nt   % going over time
       lambda_1=vp_1*dt*j;  % displacement for external protocal with vp_1
       lambda_2=vp_2*dt*j;  % displacement for external protocal with vp_2
       rold_1(N+1,1)=lambda_1;  % setting lambda(j)=lambda_1 for simulating the process _1
       rold_2(N+1,1)=lambda_2;  % setting lambda(j)=lanbda_2 for calculating _bias
       dnoise=sqrt(2*eta*kBT*dt)*randn(N,1);
       DV_1=DV1(rold_1);
       force=DV_1(1:N);         % potential gradient for simulating the process _1
       force_end_1=DV_1(N+1);   % external force for simulating the process _1
       DV_2=DV2(rold_2);        % potential gradient for simulating the process _bias
       force_end_2=DV_2(N+1);   % external force for simulating the process _bias

       DV = DV_1 - DV_2;
       DVb = DV(1:N);           % bias potential V_bias
       
       %Euler-Maruyama
       rnew(1:N,1) = rold_1(1:N) - force*dt/eta + dnoise/eta;
       W_1(j+1) = W_1(j) + force_end_1*vp_1*dt;     % work for process 1
       W_2(j+1) = W_2(j) + force_end_2*vp_2*dt;     % work for calculating _bias
       
       factor(r,j+1) = factor(r,j) + DVb'*(DVb*dt-2*dnoise)/(4*eta);	% I_bias
       P=exp(-factor(r,j+1)/kBT);
       Pfactor(r,j+1)=P;            % P_bias = exp(-beta * I_bias)
       normalization(j) = normalization(j) + P;     % summation of P_bias over realization

       % for average work
       Wave_1(j)=Wave_1(j) + W_1(j+1);
       Wave_bias(j) = Wave_bias(j) + W_2(j+1)*P;
        
       % for average displacement
       r_1(:,j)=r_1(:,j)+rnew(1:N);
       r_bias(:,j)=r_bias(:,j)+rnew(1:N)*P;
       
       % for average force
       Force_1(j) = Force_1(j) + force_end_1;
       Force_bias(j) = Force_bias(j) + force_end_2*P;
       
       %update
       rold_1 = rnew;
       rold_2 = rnew;
   end
end

%%

% Calculating the ensemble average of each observables for all time
% by /NR for _1 & _2
% and /normalization for _bias

r_1=r_1/NR;
r_2=r_2/NR;
r_bias=r_bias./normalization;

Force_1=Force_1/NR;
Force_2=Force_2/NR;
Force_bias=Force_bias./normalization;

Wave_1=Wave_1/NR;
Wave_2=Wave_2/NR;
Wave_bias=Wave_bias./normalization;


% means and variances for I_bias (factor) and P_bias (Pfactor)
mean_factor = mean(factor(:,2:end),1);
mean_Pfactor = mean(Pfactor(:,2:end),1);
var_factor = var(factor(:,2:end),1);
var_Pfactor = var(Pfactor(:,2:end),1);



%% Comparing the ensemble averages of observables with both calculation strategies and analytical result

time= [dt:dt:T];    % all timesteps for Langevin


figure;
plot(time,r_2,'b','LineWidth',2,'DisplayName','$<x>$, $V_2$');
hold on
plot(time,r_1,'Color',[0.9290 0.6940 0.1250],'LineWidth',2,'DisplayName','$<x>$, $V_1$');
hold on
plot(time,r_bias,'r:','LineWidth',2,'DisplayName','$<r>$,  $V_1$ with bias');
hold off
xlabel('time');
ylabel('average position');
legend('Location','NorthWest','Interpreter','latex');

figure;
plot(time,Force_2,'b','LineWidth',2,'DisplayName','$<x>$, $V_2$');
hold on
plot(time,Force_bias,'r:','LineWidth',2,'DisplayName','$<r>$,  $V_1$ with bias');
hold on
plot(time,Force_1,'Color',[0.9290 0.6940 0.1250],'LineWidth',2,'DisplayName','$<x>$, $V_1$');
xlabel('time');
ylabel('average external force');
legend('<F> V_2','<F> V_1 with bias','<F> V_1');
legend('Location','NorthWest','Interpreter','latex');

figure;
plot(time,Wave_2,'b','LineWidth',2,'DisplayName','$<x>$, $V_2$');
hold on
plot(time,Wave_bias,'r:','LineWidth',2,'DisplayName','$<r>$,  $V_1$ with bias');
hold on
plot(time,Wave_1,'Color',[0.9290 0.6940 0.1250],'LineWidth',2,'DisplayName','$<x>$, $V_1$');xlabel('time');
xlabel('time');
ylabel('work average');
legend('<W> V_2','<W> V_1 with bias','<W> V_1');
legend('Location','NorthWest','Interpreter','latex');



%% Show the statistics of I_bias (factor) and P_bias (Pfactor) v.s. t


figure
plot(time,mean_factor,'LineWidth',3)
set(gca,'FontSize',20);
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$\left< I_b \right>$','Interpreter','latex','FontSize',20)


figure
plot(time,mean_factor,'LineWidth',3)
set(gca,'FontSize',20,'xscale','log','yscale','log');
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$\left< I_b \right>$','Interpreter','latex','FontSize',20)


figure
plot(time,var_factor,'LineWidth',3)
set(gca,'FontSize',20);
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$Var(I_b)$','Interpreter','latex','FontSize',20)


figure
plot(time,var_factor,'LineWidth',3)
set(gca,'FontSize',20,'xscale','log','yscale','log');
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$Var(I_b)$','Interpreter','latex','FontSize',20)


figure
plot(time,mean_Pfactor,'LineWidth',3)
set(gca,'FontSize',20);
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$\left< P_b \right>$','Interpreter','latex','FontSize',20)


figure
plot(time,mean_Pfactor,'LineWidth',3)
set(gca,'FontSize',20,'xscale','log','yscale','log');
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$\left< P_b \right>$','Interpreter','latex','FontSize',20)


figure
plot(time,var_Pfactor,'LineWidth',3)
set(gca,'FontSize',20);
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$Var(P_b)$','Interpreter','latex','FontSize',20)


figure
plot(time,var_Pfactor,'LineWidth',3)
set(gca,'FontSize',20,'xscale','log','yscale','log');
xlabel('$t$','Interpreter','latex','FontSize',20)
ylabel('$Var(P_b)$','Interpreter','latex','FontSize',20)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save data

mkdir('Result')
FileName = ['_N' num2str(N,'%g') '_R' num2str(NR,'%g') '_k' num2str(k2_1,'%g') '_' num2str(k2_2,'%g') ...
    '_4k' num2str(k4_1,'%g') '_' num2str(k4_2,'%g') '_eta' num2str(eta,'%g') '_kT' num2str(kBT,'%g')...
    '_v' num2str(vp_1,'%g') '_' num2str(vp_2,'%g') ...
    '_dt' num2str(dt,'%g')];

save(['Result/DataLang' FileName '.mat'], 'N', 'eta', 'kBT', 'vp_1', 'vp_2', 'NR', 'T', 'dt', 'Nt', 'type')
if  strcmp(type,'quartic')
    save(['Result/DataLang' FileName '.mat'], 'k2_1', 'k4_1', 'k2_2', 'k4_2', '-append')
else
    disp('Wrong potential type');
    stop
end

save(['Result/DataLang' FileName '.mat'], 'time', 'r_1', 'r_2', 'r_bias', 'Wave_1', 'Wave_2', 'Wave_bias', ...
    'Force_1', 'Force_2', 'Force_bias', 'normalization', ...
    'mean_factor', 'mean_Pfactor', 'var_factor', 'var_Pfactor', '-append')



if PfactorOutput == 1       % save P_bias for all realization if PfactorOutput == 1
    index_Pb = 0:20:length(time);       % index for the coarser time grid
    index_Pb(1) = [];
    time_Pb = time(index_Pb);           % a coarser time grid for outputing Pfactor
    Pfactor_out = Pfactor(:,index_Pb);  % the output Pfactor with coarser time grid
    save(['Result/DataPb' FileName '.mat'], 'Pfactor_out', 'time_Pb','index_Pb')
end


% Save figures

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
mkdir('Result')
mkdir('Result/figures')


savefig(FigList,['Result/figures/FigLang' FileName '.fig'])


%% Plot and save the evolutions of the histogram of I_bias and P_bias

mkdir('Result/movies')
vedio = VideoWriter(['Result/movies/hist_Ib' FileName], 'MPEG-4'); % video name
vedio.FrameRate = 5; % frame speed
open(vedio)

fig = figure;

for it = Nt/1e3:Nt/1e3:Nt
    
    histogram(factor(:,it+1)/normalization(it));
    hold off
    set(gca,'FontSize',20);
    title(['time = ',num2str(time(it))])
    xlabel('$I_b$','Interpreter','latex','FontSize',20)
    ylabel('','Interpreter','latex','FontSize',20)
    frame = getframe(fig);      
    writeVideo(vedio, frame);
end
close(vedio);


vedio = VideoWriter(['Result/movies/hist_Pb' FileName], 'MPEG-4'); % video name
vedio.FrameRate = 5; % frame speed
open(vedio)

fig = figure;

for it = Nt/1e3:Nt/1e3:Nt
    
    histogram(Pfactor(:,it+1)/normalization(it));
    hold off
    set(gca,'FontSize',20);
    title(['time = ',num2str(time(it))])
    xlabel('$e^{-\beta I_b} / \mathcal{N}$','Interpreter','latex','FontSize',20)
    ylabel('','Interpreter','latex','FontSize',20)
    frame = getframe(fig);      
    writeVideo(vedio, frame);
end
close(vedio);

TClock_f = clock;
calTime_Lang = etime(TClock_f,TClock_0)     % total calculation time

save(['Result/DataLang' FileName '.mat'], 'calTime_Lang', '-append')
