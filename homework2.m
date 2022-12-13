%% Numerical Methods Homework #2
% Author: Marco Marrufo

%% Eq. (3)
clear, clc

% Simulated beta_cr values from homework #1
beta1 = 4.67;
beta2 = 7.8449;
beta3 = 11.0123;


% Circuit Parameters
f_H = 1e9;
tau = 1/(2*pi*f_H);
phi = -pi/4;

% Beta test cases -- uncomment the one you're using.
% beta = beta1/2;
% beta = (beta1+beta2)/2;
beta = (beta2+beta3)/2;


% Differential equation (3)
dxdt = @(t, x) -1/tau*x + beta/tau*cos(x+phi)^2;

% Step size
h = 1/1e12;

% Initial conditions
t0 = 0;
x0 = 0.2;


% Max simulation time & number of iterations for initialization.
tmax = 1.5e-9;
iters = int64((tmax - t0)/h)+1;

% Euler's Method

disp("Running Euler's Method")

% Starting conditions
t = t0;
x = x0;


% Initializing empty vectors of zeros.
euler_tvals = zeros(1, iters);
euler_xvals = zeros(1, iters);

i = 1;

% First step of eulers
euler_tvals(i) = t0;
euler_xvals(i) = x0;

while t < tmax
    i = i + 1;
    
    % Store previous values
    t_prev = t;
    x_prev = x;
    
    % Find next values
    t = t_prev + h;
    x = x_prev + h*dxdt(t_prev, x_prev);
    
    % Save inside vector
    euler_tvals(i) = t;
    euler_xvals(i) = x;
    
end


% RK4

disp("Running RK4")

t = t0;
x = x0;

rk4_tvals = zeros(1, iters);
rk4_xvals = zeros(1, iters);

i = 1;

rk4_tvals(i) = t0;
rk4_xvals(i) = x0;

while t < tmax
    i = i + 1;
    
    t_prev = t;
    x_prev = x;
    
    % RK4 Calculations
    k1 = dxdt(t_prev, x_prev);
    k2 = dxdt(t_prev + h/2, x_prev + h/2*k1);
    k3 = dxdt(t_prev + h/2, x_prev + h/2*k2);
    k4 = dxdt(t_prev + h, x_prev + h*k3);
    
    t = t_prev + h;
    x = x_prev + h/6*(k1+2*k2+2*k3+k4);
    
    rk4_tvals(i) = t;
    rk4_xvals(i) = x;
end

%
figure('Position', [650 350 900 600])
plot(euler_tvals, euler_xvals, 'k');
hold on
plot(rk4_tvals, rk4_xvals);
title("Numerical Methods for $$\frac{dx}{dt} = -\frac{1}{\tau}x" ...
    + "+ \frac{\beta}{\tau}\cos^2{[x+\phi]}$$" + ", $$\beta$$ = " + beta, 'interpreter','latex');
legend(["Euler's", "RK4"], "Location", "best")
xlim([t0-h, tmax+h])
xlabel("t")
ylabel("x")


%% Eqs. (8) and Eq (9)
clear, clc

% Circuit Parameters
f_H = 1e9;
f_L = 1e6;
phi = -pi/4;

tau = 1/(2*pi*f_H);
theta = 1/(2*pi*f_L);

beta_cr = -1/sin(2*phi);

% Different test cases for beta -- change gain as needed
% (1) gain = 0.99
% (2) gain = 1.01
% (3) gain = 3
gain = 3;
beta = gain*beta_cr;



xdot = @(x, y) -1/tau*x - 1/tau*y + beta/tau*cos(x + phi)^2;
ydot = @(x) 1/theta*x;

% Step size
h = 1/1e11;

% Initial conditions
t0 = 0;
x0 = 0;
y0 = 0;


% t max
tmax = 1.1e-5;
iters = int64((tmax - t0)/h)+1;

% Euler's Method

disp("Running Euler's Method")

t = t0;
x = x0;
y = y0;


euler_tvals = zeros(1, iters);
euler_xvals = zeros(1, iters);
euler_yvals = zeros(1, iters);

i = 1;

euler_tvals(i) = t0;
euler_xvals(i) = x0;
euler_yvals(i) = y0;

while t < tmax

    i = i + 1;
    
    t_prev = t;
    x_prev = x;
    y_prev = y;
    
    t = t_prev + h;
    x = x_prev + h*xdot(x_prev, y_prev);
    y = y_prev + h*ydot(x_prev);
    
    euler_tvals(i) = t;
    euler_xvals(i) = x;
    euler_yvals(i) = y;
    
end



% RK4

disp("Running RK4")

% Initial conditions
t = t0;
x = x0;
y = y0;

% Initializing empty vectors
rk4_tvals = zeros(1, iters);
rk4_xvals = zeros(1, iters);
rk4_yvals = zeros(1, iters);

i = 1;

rk4_tvals(i) = t0;
rk4_xvals(i) = x0;
rk4_yvals(i) = y0;

while t < tmax
    i = i + 1;
    
    t_prev = t;
    x_prev = x;
    y_prev = y;
    
    % RK4 Calculations
    k1x = xdot(x_prev, y_prev);
    k1y = ydot(x_prev);
    
    k2x = xdot(x_prev + h/2*k1x, y_prev + h/2*k1y);
    k2y = ydot(x_prev + h/2*k1x);
    
    k3x = xdot(x_prev + h/2*k2x, y_prev + h/2*k2y);
    k3y = ydot(x_prev + h/2*k2x);

    k4x = xdot(x_prev + h*k3x, y_prev + h*k3y);
    k4y = ydot(x_prev + h*k3x);
    
    t = t_prev + h;
    x = x_prev + h/6*(k1x+2*k2x+2*k3x+k4x);
    y = y_prev + h/6*(k1y+2*k2y+2*k3y+k4y);
    
    rk4_tvals(i) = t;
    rk4_xvals(i) = x;
    rk4_yvals(i) = y;
end

% Plots

mesh_x = linspace(min(rk4_xvals), max(rk4_xvals), 20);
mesh_y = linspace(min(rk4_yvals), max(rk4_yvals), 20);

figure('Position', [650 350 900 600])
[X, Y] = meshgrid(mesh_x, mesh_y);

u = zeros(size(X));
v = zeros(size(X));

for i = 1:numel(X)
   u(i) = xdot(X(i), Y(i));
   v(i) = ydot(X(i));
end

h = quiver(X, Y, u, v, 'r'); figure(gcf)
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

hold on

titletext = strcat("Numerical Method for Eq. (8) and Eq. (9), $$\beta =  $$", " ", num2str(gain),  "$$\beta_{Cr}$$");
title(titletext ,'interpreter','latex')
plot(euler_xvals, euler_yvals, 'k', 'LineWidth', 0.5);
plot(rk4_xvals, rk4_yvals, 'LineWidth', 1.5);
plot(rk4_xvals(1),rk4_yvals(1),'mo','LineWidth', 1.5) % starting point
plot(rk4_xvals(end),rk4_yvals(end),'g*','LineWidth', 1.5) % ending point
legend(["Euler's", "RK4", "Start", "End"], "Location", "best")
xlabel("x")
ylabel("y")
hold off

figure()
plot(rk4_tvals, rk4_xvals)

%% Eqs. (13) and (14)
clear, clc

% Circuit Parameters
omega_0 = 1e9;
delta_omega = 100e6;
phi = -pi/4;



beta_cr = -1/sin(2*phi);

% Different test cases for beta -- change gain as needed
% (1) gain = 0.99
% (2) gain = 1.01
% (3) gain = 3
gain = 0.99;

beta = gain*beta_cr;



xdot = @(x, y) delta_omega*(-x - y + beta*cos(x + phi)^2);
ydot = @(x) omega_0^2/delta_omega*x;

% Step size
h = 1/1e12;

% Initial conditions
t0 = 0;
x0 = 0.2;
y0 = 0.2;


% t max
tmax = 1.1e-5;
iters = int64((tmax - t0)/h)+1;

% Euler's Method

disp("Running Euler's Method")

t = t0;
x = x0;
y = y0;


euler_tvals = zeros(1, iters);
euler_xvals = zeros(1, iters);
euler_yvals = zeros(1, iters);

i = 1;

euler_tvals(i) = t0;
euler_xvals(i) = x0;
euler_yvals(i) = y0;

while t < tmax

    i = i + 1;
    
    t_prev = t;
    x_prev = x;
    y_prev = y;
    
    t = t_prev + h;
    x = x_prev + h*xdot(x_prev, y_prev);
    y = y_prev + h*ydot(x_prev);
    
    euler_tvals(i) = t;
    euler_xvals(i) = x;
    euler_yvals(i) = y;
    
end

% RK4

disp("Running RK4")

% Initial conditions
t = t0;
x = x0;
y = y0;

% Initializing empty vectors
rk4_tvals = zeros(1, iters);
rk4_xvals = zeros(1, iters);
rk4_yvals = zeros(1, iters);

i = 1;

rk4_tvals(i) = t0;
rk4_xvals(i) = x0;
rk4_yvals(i) = y0;

while t < tmax
    i = i + 1;
    
    t_prev = t;
    x_prev = x;
    y_prev = y;
    
    % RK4 Calculations
    k1x = xdot(x_prev, y_prev);
    k1y = ydot(x_prev);
    
    k2x = xdot(x_prev + h/2*k1x, y_prev + h/2*k1y);
    k2y = ydot(x_prev + h/2*k1x);
    
    k3x = xdot(x_prev + h/2*k2x, y_prev + h/2*k2y);
    k3y = ydot(x_prev + h/2*k2x);

    k4x = xdot(x_prev + h*k3x, y_prev + h*k3y);
    k4y = ydot(x_prev + h*k3x);
    
    t = t_prev + h;
    x = x_prev + h/6*(k1x+2*k2x+2*k3x+k4x);
    y = y_prev + h/6*(k1y+2*k2y+2*k3y+k4y);
    
    rk4_tvals(i) = t;
    rk4_xvals(i) = x;
    rk4_yvals(i) = y;
end

% Plots

mesh_x = linspace(min(rk4_xvals), max(rk4_xvals), 20);
mesh_y = linspace(min(rk4_yvals), max(rk4_yvals), 20);

%figure('Position', [650 350 900 600])
[X, Y] = meshgrid(mesh_x, mesh_y);

u = zeros(size(X));
v = zeros(size(X));

for i = 1:numel(X)
   u(i) = xdot(X(i), Y(i));
   v(i) = ydot(X(i));
end

%h = quiver(X, Y, u, v, 'r'); figure(gcf)
%h.Annotation.LegendInformation.IconDisplayStyle = 'off';

%hold on

%titletext = strcat("Numerical Method for Eq. (8) and Eq. (9), $$\beta =  $$", " ", num2str(gain),  "$$\beta_{Cr}$$");
%title(titletext ,'interpreter','latex')
% plot(euler_xvals, euler_yvals, 'k', 'LineWidth', 0.5);
% comet(rk4_xvals, rk4_yvals, 'LineWidth', 1.5);
comet3(rk4_tvals, rk4_xvals, rk4_yvals);
%plot(rk4_xvals(1),rk4_yvals(1),'mo','LineWidth', 1.5) % starting point
%plot(rk4_xvals(end),rk4_yvals(end),'g*','LineWidth', 1.5) % ending point
%legend(["Euler's", "RK4", "Start", "End"], "Location", "best")
%xlabel("x")
%ylabel("y")
%hold off

%figure()
%plot(rk4_tvals, rk4_xvals)