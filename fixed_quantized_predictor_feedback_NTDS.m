%Reset
close all
clear all 
clc
global  U L k

t_y= 0; % initial time t_0
Ts=10;% time limit 

L=1;
k=1;
D=1; % The delay constant 
Nx= 100; % the number of points of space discretization (number of mesh points)
x= linspace(0,D,Nx); % the space discretization
dx= x(2) - x(1); % the space discretization step

% Time discretization :
t_f= t_y+Ts; % final time for integration on intervale [t_i,t_i+Ts]
Nt= 100; % the number of points of time discretization (time)
howfar= t_f/Nt;
t= t_y; % the initialization time

coef= 0.65; % coef depends on lambda_max to satisfy CFL condition
timestep= coef*dx;% Time discretization's step

% Initialization of the solutions by zeros
Yevol1= zeros(Nt,Nx); % X(t,x)

tout_y= zeros(Nt,1); %initialize the time

Y= InitialConditions(x); %  Y(t0,x)
%Y= zeros(3,Nx);  
Yevol1(1,:)= Y(1,:); %  X_1(0), ..., X_1(0) Nx times
Yevol3(1,:)= Y(3,:); %  u(0,x)


mu_t=100;%mu_t=0.1
M=2;
Delta=M/100;
Mbar1= 2;
t0=.0;
U_values = zeros(Nt,1); % Initialize array to store U(t)
%mu_values = zeros(Nt,1); % Initialize array to store mu(t)

tic
% Defining the structure of the original system
sol= setup(1,@DefineTransportEquation,t,x,Y,'LxF',[],@TransportBoundaryConditions); 

Integral = zeros(size(x));
if  t_y<=t0
            U = 0;
end % 
for m= 2:Nt  
    sol= hpde(sol,howfar,timestep); % The solver
    t_y= sol.t; % The current time 
    Y= sol.u; % The solutions of the original system   

 if  t_y<=t0
            U = 0;
 else 
      % Compute P_mu(x) using the predictor equation
    p_mu(1) = mu_t * quantizer(Y(1,1), mu_t, M, Delta); % Boundary condition for p_mu at x = 0
    for i = 1:Nx-1
         f_p = sign(p_mu(i)) * (p_mu(i)^2 / sqrt(1 + p_mu(i)^2)) + mu_t * quantizer(Y(3,i),mu_t, M, Delta); 
        p_mu(i+1) = p_mu(i) + dx * f_p; 
    end
    P_mu = p_mu(end); % P_mu = p_mu(D)
    U = -(k + L) * P_mu; % Feedback control
 end

    tout_y(m)= t_y; % Storage of the solutions to use later to plot the solutions    
    Yevol1(m,1)= max(Y(1,1),0); % X 
    Yevol3(m,:)= Y(3,:); % u(t,x)   
    U_values(m) = U;
end  
toc
% Initialize storage for norms of X and u
norm_X = zeros(Nt, 1);
sup_norm_u = zeros(Nt, 1);
for m = 1:Nt
    % Compute norms of X and sup norm of u at each time step
    norm_X(m) = abs(Yevol1(m,1));
    sup_norm_u(m) = max(abs(Yevol3(m,:)));
end

% Plot state evolution of ODE over time
Nfig= 0;    
%Figure 1 : The solution X(t)
Nfig= Nfig+1; figure(Nfig)
% Plot X(t)
subplot(2,2,[1,2]);
ax=gca;
ax.FontSize=35;
plot(tout_y, Yevol1(:,1),'k','LineWidth', 5 ); % plot X_1(t)
xlabel('$t$', 'Interpreter', 'latex','LineWidth', 5, 'FontSize', 30);
ylabel('$X(t)$', 'Interpreter', 'latex','LineWidth', 5, 'FontSize', 30);
set(gca, 'FontSize', 35);
hold on 
grid on 
% Modify the grid properties
ax = gca; % Get the current axes
ax.GridLineStyle = '-'; % Style of the grid lines ('-' for solid lines, '--' for dashed, etc.)
ax.GridColor = [0, 0, 0]; % Color of the grid lines (black in this example)
ax.GridAlpha = 0; % Transparency of the grid lines (1 = opaque, 0 = transparent)
ax.LineWidth = 3; % Thickness of the grid lines

% Plot coupled solution of transport equation and ODE
subplot(2,2,[3,4]);
pstep = 3; % Increase the spacing between spatial grid points in the plot
tstep = 3; % Increase the spacing between time steps in the plot
[Y, Z] = meshgrid(x(1:pstep:end), howfar*(0:tstep:Nt-1));
U_plot = Yevol3(1:pstep:end, 1:tstep:end);
mesh(Y, Z, U_plot, 'LineWidth', 1,'edgecolor', 'black');
view(83,10);
ax=gca;
ax.FontSize=35;
xlabel('$x$', 'interpreter', 'latex', 'LineWidth', 5,'FontSize', 30);
ylabel('$t$', 'interpreter', 'latex', 'LineWidth', 5,'FontSize', 30);
zlabel('$u(x,t)$', 'interpreter', 'latex', 'LineWidth', 5,'FontSize', 30);
set(gca, 'FontSize', 35);
hold on 
grid on 
% Modify the grid properties
ax = gca; % Get the current axes
ax.GridLineStyle = '-'; % Style of the grid lines ('-' for solid lines, '--' for dashed, etc.)
ax.GridColor = [0, 0, 0]; % Color of the grid lines (black in this example)
ax.GridAlpha = 0; % Transparency of the grid lines (1 = opaque, 0 = transparent)
ax.LineWidth = 3; % Thickness of the grid lines


%%%%%%%%%%%%%%%%%%%%%%%% Figure 4
Nfig = Nfig + 1;
figure(Nfig);
ax=gca;
ax.FontSize=25;
plot(tout_y, sup_norm_u+norm_X,'k','LineWidth', 5);
xlabel('$t$','Interpreter','latex','FontSize',30);
ylabel('$|X(t)|+\|u(t)\|_{\infty}$','Interpreter','latex','FontSize',30);
set(gca, 'FontSize', 35);
grid on 
hold on 
% Modify the grid properties
ax = gca; % Get the current axes
ax.GridLineStyle = '-'; % Style of the grid lines ('-' for solid lines, '--' for dashed, etc.)
ax.GridColor = [0, 0, 0]; % Color of the grid lines (black in this example)
ax.GridAlpha = 0; % Transparency of the grid lines (1 = opaque, 0 = transparent)
ax.LineWidth = 3; % Thickness of the grid lines

% Definition of the inital condition of ODE-PDE system
function Y = InitialConditions(x)
      Y(1,:)= 0*x + 1*1; %X_1(0)
      Y(3,:)= 0; % u(0,x)%
end   
% Definition of the linear ODE-PDE system
function Y_t= DefineTransportEquation(t,x,Y,Y_x) 
    Y_t(1,:)= sign(Y(1,:)) * (Y(1,:).^2 / sqrt(1 + Y(1,:).^2)) + Y(3,1); % ODE X(t)
    
    Y_t(3,:)=  Y_x(3,:)+0*Y(3,:); % PDE u(t,x)
end
% Definition of the boundary condition of the ODE-PDE system
function [YL,YR]= TransportBoundaryConditions(t,YLex,YRex) 
    global U
    YL(1)= YLex(1); % No left boundary X_1(t)
    YR(1)= YRex(1); % No left boundary X_1(t)
    YL(3)= YLex(3); % No left boundary u(t,x)
    YR(3)= U;  % Control
end 
% quantizer
function quantized_u = quantizer(u,mu_t, M, Delta)
    if u/mu_t >= M
        quantized_u= M;
    elseif u/mu_t <= -M
        quantized_u= -M;
    else
        quantized_u= Delta* floor(u/(Delta*mu_t) + 0.5);
    end
end