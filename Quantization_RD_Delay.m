% Heat equation  in 1D 
clc
clear all
close all
format long

D = 1; % The delay
lambda = 12; % reaction coefficient
%%%%%%%%%%%%%%% time discretization   %%%%%%%%%%%%%%%%%%%
tf = 10; 
nt = 1000;
dt = tf/nt;
t = 0:dt:tf;
%%%%%%%%%%%%%% space discretization %%%%%%%%%%%%%%%%%%%%
a = 0;
b = 1;
nx = 100;%
dx = b/nx;
x = a:dx:b;
y= a:dx:b;
%-----------------------------------------------------------%
zL = @(t) 0;% left boundary
zR = @(t) 0;% right boundary
uinitial = @(x) u0(x);% initial condition
vinitial = @(x) x-x+5;% initial condition
%%%%%%%%%%%%%% initialization of u(t,x) and v(t,x) %%%%%%%%%%%%%%%%%%%%
usol(1,:) = uinitial(x); % Initialize u(t,x)
vsol(1,:) = vinitial(x); % Initialize v(t,x)
vsol(1,1) = usol(1,end); % Enforce boundary condition u(t,1) = v(t,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coef1 = dt/(dx)^2;
coef2 = dt/(D*dx);
A = zeros(nx+1+nx+1);
A(1,1) = 1;
for i = 2:nx        
    A(i,i-1:i+1) = [-coef1 1+2*coef1-lambda*dt -coef1];
end
A(nx+1,nx+1:nx+2) = [-1 1];
A(nx+2,nx+2:nx+3) = [(1+coef2) -coef2];
for i = nx+3:nx+1+nx        
    A(i,i-1:i+1) = [0.5*coef2 1 -0.5*coef2];
end
A(end,end) = 1;

B = zeros(nx+1+nx+1,1);
B(end,1) = 1;
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%% Calculate the kernel k(x,y) %%%%%%%%%%%%%%%%%%%%
for i=1:nx+1
    for j=1:i
        kernel_k_xy(i,j)= -lambda*x(j)*besseli(1,sqrt(lambda*(x(i)^2-x(j)^2)))/(sqrt(lambda*(x(i)^2-x(j)^2)));
    end
    kernel_k_xy(i,i) = -lambda/2*x(i);
end
%%%%%%%%%%%%%%%% Calculate the kernels gamma_y(x-y,1) %%%%%%%%%%%%%%%%%%%%
n_inf = 30; % max=162 the number of terms in the serie 
Gamma = zeros(nx+1,nx+1);
Gamma_y = zeros(nx+1,nx+1);
for n=1:n_inf  
    int_x(1) = 0;
    for i=2:nx+1
        int_x(i) = trapz(x(1:i),kernel_k_xy(i,1:i).*sin(n*pi*x(1:i)));
    end
    for i=1:nx+1
        for j=1:i
            Gamma_y(i,j) = Gamma_y(i,j) + 2*pi*n*(-1)^n*int_x(end)*exp(D*(lambda-pi^2*n^2)*(x(i)-y(j))); % gamma_y(x-y,1)
        end
    end
end
%%%%%%%%%%%%%%%% Calculate gamma(x,y) %%%%%%%%%%%%%%%%%%%%%%%%%
n_inf = 30; % max=162 the number of terms in the serie 
for n=1:n_inf  
    int_x(1) = 0;
    for i=2:nx+1
        int_x(i) = trapz(x(1:i),kernel_k_xy(i,1:i).*sin(n*pi*x(1:i)));
    end
    for i=1:nx+1    
        for j=1:nx+1
            Gamma(i,j) = Gamma(i,j) + 2*int_x(end)*exp(D*(lambda-pi^2*n^2)*x(i)).*sin(n*pi*y(j)); % gamma(x,y)   
        end
    end 
end

% %%%%%%%%%%%%%%%%%%%%  Control Backstepping initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % control = -trapz(x,Gamma(end,:).*usol(1,:)) - D*trapz(y,Gamma_y(end,:).*vsol(1,:));
% vsol(1,end) = control;% v(t_0,1) = U(t)
%%%%%%%%%%%%%% Initialization of the the norm vectors %%%%%%%%%%%%%%%%%%%%%%%%
L2normusol = L2norm_function(x,usol(1,:)); % ||u(t_0,.)||_L^2
L2normvsol = L2norm_function(y,vsol(1,:)); % ||v(t_0,.)||_L^2

% Parameters for quantized control
M = 1; 
Delta =M/100; 
Mbar = .1; 
mu0 = .1;
Omega = 0.3;
t0 = 0; 
T = 1; 
tau = 1;
sigma1=lambda-pi.^2;
M3=sqrt(L2norm_function(x,Gamma(end,:)))+max(abs(Gamma_y(end,:)));

mu_values = zeros(nt,1); % Initialize array to store mu(t)
mu_values(1)=NaN;

 N_max = 30; % Large value for convergence
 sum_term = 0; % Initialize sum
for n = 1:N_max  % Start from n = 1 (since n=0 term is 0)
    sum_term = sum_term + (n^2 * pi^2) / ((lambda - n^2 * pi^2)^2);
end
% Compute G
G = 4 * sqrt(sum_term);
Mbar1 = max(G+1,sqrt(2));
Unom_values = zeros(size(usol));
for i = 2:nt+1 % time loop 
%%%%%% solution %%%%%%%%%% 
    Z= -[usol(i-1,:) vsol(i-1,:)]';
    Z(1) = 0;
    Z(nx+1) = 0;
    Z(end) = 0;
    % -Z^n +As - B Control=0
    current_time = t(i);
    mu_t = mu(current_time, t0, tau, T, Mbar1, mu0, Omega, sigma1);
    mu_values(i) = mu_t; % Store mu_t in the correct index
     if  current_time<=t0
        U=@(s) 0;
     else 
         Unom=@(s) trapz(x,Gamma(end,:).*s(1:nx+1,1)') - D*trapz(y,Gamma_y(end,:).*s(nx+2:end,1)');
    % Unom=@(s) trapz(x,Gamma(end,:)*mu_t.*quantizer(s(1:nx+1,1)', mu_t, M, Delta)) - D*trapz(y,Gamma_y(end,:)*mu_t.*quantizer(s(nx+2:end,1)', mu_t, M, Delta));%;trapz(x,Gamma(end,:).*s(1:nx+1,1)') - D*trapz(y,Gamma_y(end,:).*s(nx+2:end,1)');
         U=@(s) mu_t*quantizer(Unom(s), mu_t, M, Delta);
     end 
 
     %trapz(x,Gamma(end,:).*mu_t*quantizer(s(1:nx+1,1)', mu_t, M, Delta)) - D*trapz(y,Gamma_y(end,:).*quantizer(s(nx+2:end,1)', mu_t, M, Delta));
    Implicit_Numerical_Scheme = @(s) A*s - B*U(s) + Z;%      (kernel_integral*(dx*kernel_term1*x(1:nx+1,1)+ dy*kernel_term2*x(nx+2:end,1)))
    options = optimset('TolFun',1e-100,'Display','off');
    sol = fsolve(Implicit_Numerical_Scheme,[usol(i-1,:) vsol(i-1,:)]',options);%
    usol(i,:) = sol(1:nx+1)';
    vsol(i,:) = sol(nx+2:end)';
    %%%%%%%%%%%%%%%% Extract and rename solutions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% calculating the the norm vectors at time t %%%%%%%%%%%%%%
    L2normusol = [L2normusol, L2norm_function(x,usol(i,:))];
    sup_norm_v = sup_norm_function(vsol(i,:));
    L2normvsol = [L2normvsol, sup_norm_v];% sup norm calculation
end
%--------------------------------------- Plots ---------------------------%
pstep = 3; % Increase the spacing between spatial grid points in the plot
tstep = 3; % Increase the spacing between time steps in the plot
%%%%%%%%%%%%%%%%%%%%%%%%%% Plot for u(t,x)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
%subplot(2,1,1);
[Y, Z] = meshgrid(x(1:pstep:end), t(1:tstep:end)); % Create meshgrid for x and t
U_plot = usol(1:tstep:end, 1:pstep:end); % Extract u(t,x) with spacing
mesh(Y, Z, U_plot, 'LineWidth', 1, 'edgecolor', 'black'); % Plot u(t,x)
view(83,10); % Set view angle
ax = gca;
ax.FontSize = 30; % Set font size
xlabel('$x$', 'interpreter', 'latex', 'LineWidth', 5, 'FontSize', 35); % Label x-axis
ylabel('$t$', 'interpreter', 'latex', 'LineWidth', 5, 'FontSize', 35); % Label y-axis
zlabel('$u(x,t)$', 'interpreter', 'latex', 'LineWidth', 5, 'FontSize', 35); % Label z-axis
set(gca, 'FontSize', 35); % Set axis font size
hold on
grid on
% Modify the grid properties
ax.GridLineStyle = '-'; % Solid grid lines
ax.GridColor = [0, 0, 0]; % Black grid lines
ax.GridAlpha = 0.5; % Semi-transparent grid lines
ax.LineWidth = 3; % Thickness of the grid lines
% Plot for v(t,y)
figure(2)
% Create meshgrid for y and t
[Y, Z] = meshgrid(y(1:pstep:end), t(1:tstep:end)); % Mesh for v(t,y)
V_plot = vsol(1:tstep:end, 1:pstep:end); % Extract v(t,y) with spacing
% Plot the mesh
mesh(Y, Z, V_plot, 'LineWidth', 1, 'edgecolor', 'black'); 
view(83,10); % Set view angle
% Axis settings
ax = gca;
ax.FontSize = 35;
xlabel('$x$', 'interpreter', 'latex', 'LineWidth', 5, 'FontSize',35);
ylabel('$t$', 'interpreter', 'latex', 'LineWidth', 5, 'FontSize', 35);
zlabel('$v(x,t)$', 'interpreter', 'latex', 'LineWidth', 5, 'FontSize', 35);
set(gca, 'FontSize', 35);
hold on;
grid on;
% Modify grid properties
ax.GridLineStyle = '-';
ax.GridColor = [0, 0, 0]; % Black grid lines
ax.GridAlpha = 0.5;
ax.LineWidth = 3;
%**Overlay vsol(:,end) in red**
t_vals = t(1:tstep:end); % Extract corresponding time values
y_last = y(end) * ones(size(t_vals)); % Keep y-position at the last column
v_last = vsol(1:tstep:end, end); % Extract last column of vsol
% Plot as a red 3D line over the mesh
plot3(y_last, t_vals, v_last, 'r', 'LineWidth', 5); 
hold off;
%%%%%%%%%%%%%% figure 3 %%%%%%%%%%%%%%%%%%
figure(3)
ax=gca;
%plot(t,sqrt(L2normusol)+L2normvsol,'k','LineWidth',3);
plot(t,sqrt(L2normusol)+L2normvsol,'k',t, (M*Mbar)*mu_values,'k--','LineWidth',5);
hold on;  
grid on
xlabel('t','interpreter','latex')
ylabel('$\|u(\cdot,t)\|_2+\|v(\cdot,t)\|_{\infty}$','interpreter','latex')%'$\frac{M_{2}}{M_3} M\mu(t)$',
set(gca, 'FontSize', 35); 
% %%%%
figure(4);
ax=gca;
ax.FontSize=30;
plot(t, vsol(:,end),'k','LineWidth',5);
xlabel('$t$','Interpreter','latex','LineWidth', 5,'FontSize',35);
ylabel('$U(t)$','Interpreter','latex','LineWidth', 5,'FontSize',35);
set(gca, 'FontSize', 35);
hold on 
grid on 
%Modify the grid properties
ax.GridLineStyle = '-'; % Solid grid lines
ax.GridColor = [0, 0, 0]; % Black grid lines
ax.GridAlpha = 0.5; % Semi-transparent grid lines
ax.LineWidth = 3; % Thickness of the grid lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L2norm = L2norm_function(x,L)
    b = size(L,2);   
    dx = x(2)-x(1);
    L2norm = 0;
    for k=1:b %  computing a sort of Riemman sum to approximate the integral
        L2norm = L2norm + L(1,k)*L(1,k);  
    end
    L2norm = L(1,1)*L(1,1) + dx*L2norm;
end
function sup_norm = sup_norm_function( v)
    % Compute the sup norm (infinity norm) of v
    sup_norm = max(abs(v)); % Maximum absolute value of v
end
function u0 = u0(x)
    u0 = zeros(size(x));
    for n = 1:3
        u0 = u0 + (sqrt(2) / n) * sin(n * pi * x);
    end
    u0 = u0 + 3 * (x.^2 - x.^3);
end
% Quantizer
function quantized_u = quantizer(u, mu_t, M, Delta)
    if u / mu_t >= M
        quantized_u = M;
    elseif u / mu_t <= -M
        quantized_u = -M;
    else
        quantized_u = Delta * floor(u / (Delta * mu_t) + 0.5);
    end
end
%Switching parameter
function mu_t = mu(t, t0, tau, T, Mbar1, mu0, Omega, sigma1)
    if t <= t0
        % Case 1: t is within the delay period
        j = 1; % Initialize the time step index
        while (j - 1) * tau > t || t > j * tau
            j = j + 1;
        end
        % Compute mu_t for the delay period
        mu_t = Mbar1 * exp(2 * sigma1 * tau * (j)) * mu0;

    elseif t0 < t && t <= t0 + T
        % Case 2: t is within the first switching period after t0
        mu_t = mu(t0, t0, tau, T, Mbar1, mu0, Omega, sigma1);

    else
        % Case 3: t is beyond t0 + T
        i = 2; % Initialize the switching period index
        while ~(t0 + (i - 1) * T < t && t <= t0 + i * T)
            i = i + 1;
        end
        % Compute mu_t for the switching region
        mu_t = Omega * mu(t0 + (i - 1) * T, t0, tau, T, Mbar1, mu0, Omega, sigma1);
    end
end