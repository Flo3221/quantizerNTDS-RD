% Reset
close all
clear all 
clc
global  L k  D U  Mbar1 Mbar mu0 Omega t0 T M Delta t_f tau

t_y= 0; % initial time t_0
Ts= 10;% time limit 

L=1;
k=1;

D=1; % The delay constant 
L=D; % the right bound of the domain  
Nx= 100; % the number of points of space discretization (number of mesh points)
x= linspace(0,L,Nx); % the space discretization
dx= x(2) - x(1); % the space discretization step

% Time discretization :
t_f= t_y+Ts; % final time for integration on intervale [t_i,t_i+Ts]
Nt= 400; % the number of points of time discretization (time)
howfar= t_f/Nt;
t= t_y; % the initialization time

coef= 0.5; %coef depends on lambda_max to satisfy CFL condition
timestep= coef*dx; %Time discretization's step

% Initialization of the solutions by zeros
Yevol1= zeros(Nt,Nx); % X_1(t,x)
tout_y= zeros(Nt,1); %initialize the time

Y= InitialConditions(x); %  Y(t0,x)
Yevol1(1,:)= Y(1,:); %  X_1(0), ..., X_1(0) Nx times
Yevol3(1,:)= Y(3,:); %  u(0,x)

M=2;
Delta=M/100;
Mbar1= 1;
Mbar=.6;
mu0 = 1;
Omega=.63;
t0=0;
T=2;
tau=1;

U_values = zeros(Nt,1); % Initialize array to store U(t)
mu_values = zeros(Nt,1); % Initialize array to store mu(t)
mu_values(1)=NaN;
qX1_values=zeros(Nt,1);


% Defining the structure of the original system
sol= setup(1,@DefineTransportEquation,t,x,Y,'LxF',[],@TransportBoundaryConditions); 
tic
Integral = zeros(size(x));
if  t_y<=t0
            U = 0;
end % 
for m= 2:Nt  
    sol= hpde(sol,howfar,timestep); %The solver
    t_y= sol.t; %The current time 
    Y= sol.u; %The solutions of the original system in the current time
    mu_t=mu(t_y);
    
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

    qX1=mu_t*quantizer(Y(1,1),mu_t, M, Delta);
    qu=mu_t*quantizer(Y(3,:),mu_t, M, Delta);  

    tout_y(m)= t_y; % Storage of the solutions to use later to plot the solutions    
    Yevol1(m,1)= max(Y(1,1),0); % X1   
    Yevol3(m,:)= Y(3,:); % u(t,x)  
    U_values(m) = U;
    mu_values(m) = mu_t; % Store mu(t)
    qX1_values(m)=qX1;
    qu_values(m,:)=qu;
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
Y1=quantizer(Yevol1(1,1),mu(t0), M, Delta);
Y3=quantizer(Yevol3(3,1),mu(t0), M, Delta);
Mm=(M*Mbar-Delta)*mu(t0);
% Plot state evolution of ODE over time
Nfig= 0;    
%Figure 1 : The solution X(t)
Nfig= Nfig+1; figure(Nfig)
% Plot X(t)
subplot(2,2,[1,2]);
plot(tout_y, Yevol1(:,1),'k','LineWidth', 5 ); % plot X_1(t)
ax=gca;
ax.FontSize=30;
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
% Overlay u(D,t) as a red 3D line
hold on;
t_vals = tout_y;                  % time values
x_vals = L * ones(size(t_vals)); % x = D = L
u_vals = Yevol3(:, end);         % u(D,t) is the last column of Yevol3
plot3(x_vals, t_vals, u_vals, 'r', 'LineWidth', 4); % red 3D line

ax=gca;
ax.FontSize=30;
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

%Figure $|X(t)|+\|u(t)\|_{\infty}$
Nfig = Nfig + 1;
figure(Nfig);
ax=gca;
ax.FontSize=30;
plot(tout_y, sup_norm_u+norm_X,'k',tout_y, M*Mbar*mu_values,'k--','LineWidth',5);
xlabel('$t$','Interpreter','latex','LineWidth', 5,'FontSize',30);
ylabel('$|X(t)|+\|u(t)\|_{\infty}$','Interpreter','latex','LineWidth', 5,'FontSize',30);
set(gca, 'FontSize', 35);
hold on 
grid on 
% Modify the grid properties
ax = gca; % Get the current axes
ax.GridLineStyle = '-'; % Style of the grid lines ('-' for solid lines, '--' for dashed, etc.)
ax.GridColor = [0, 0, 0]; % Color of the grid lines (black in this example)
ax.GridAlpha = 0; % Transparency of the grid lines (1 = opaque, 0 = transparent)
ax.LineWidth = 3; % Thickness of the grid lines
% %%%%%%%%%%%%%%%%%%%%%%%% Figure 2 : The solution X(t)
% Nfig= Nfig+1; figure(Nfig)
% % Plot X(t)
% subplot(2,2,[1,2]);
% plot(tout_y, max(qX1_values(:,1),0),'k','LineWidth', 3 ); % plot X_1(t)
% ax=gca;
% ax.FontSize=30;
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 30);
% ylabel('$\mu(t)q( \frac{X(t)}{\mu(t)})$', 'Interpreter', 'latex', 'FontSize', 30);
% set(gca, 'FontSize', 35); 
% hold on 
% grid on 
% % Plot coupled solution of transport equation and ODE
% subplot(2,2,[3,4]);
% pstep = 3; % Increase the spacing between spatial grid points in the plot
% tstep = 3; % Increase the spacing between time steps in the plot
% [Y, Z] = meshgrid(x(1:pstep:end), howfar*(0:tstep:Nt-1));
% U_plot = qu_values(1:pstep:end, 1:tstep:end);
% mesh(Y, Z, U_plot, 'LineWidth', 1,'edgecolor', 'black');
% view(83,10);
% ax=gca;
% ax.FontSize=30;
% xlabel('$x$', 'interpreter', 'latex', 'FontSize', 30);
% ylabel('$t$', 'interpreter', 'latex', 'FontSize', 30);
% zlabel('$\mu(t)q( \frac{u(x,t)}{\mu(t)})$', 'interpreter', 'latex', 'FontSize', 30);
% set(gca, 'FontSize', 35); % Ajustez la taille de la police des axes ici
% hold on 
% grid on 
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
%Switching parameter
function mu_t = mu(t)
global  L  Mbar1 mu0 Omega t0 T t_f tau; 
       if t<=t0 
          j=1;
          while (j-1)*tau > t && t>j*tau && j<=floor(t0/tau)+1
          j=j+1;
          end 
          mu_t=Mbar1*exp(2 * L * tau * j) * mu0;
      elseif t0<t && t<= t0+T           
          mu_t=mu(t0);
      else        
           i=2;
            test=t0+(i-1)*T<t && t<=t0+i*T;
          while test==0 && i<=floor((t_f-t0)/T)
          i=i+1;
          test=t0+(i-1)*T<t && t<=t0+i*T;
          end
          mu_t=Omega*mu(t0+(i-1)*T);
      end
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