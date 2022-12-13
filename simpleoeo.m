%% RK4 for DDE Modeling of a Simple OEO w 4km Delay Line
% Author: Marco Marrufo
% Code based off of DDE numerical integration slides by Professor Chembo
clear, clc

% System parameters
fh = 6.4e6;
fl = 531.5e3;

tau = 1/(2*pi*fh);
theta = 1/(2*pi*fl);

n = 1.5;
L = 40;
c = 3e8;
Td = n*L/c;


% Tunable parameters
beta2 = (1+(Td+tau)/theta)^(1/2);
beta1 = (1+(Td/(2*pi*theta)-2*pi*tau/Td)^2)^(1/2);
beta0 = 0.99;

disp("\beta_2 is: " + beta2)
disp("\beta_1 is: " + beta1)
disp("\beta_0 is: " + beta0)

beta = beta2;
alpha = -0.5;


% Flow Equations
D = @(x) ((x).*(x > 0));
xdot = @(x, xt, y) (-1/tau*x - 1/(tau*theta)*y + beta/tau*D(xt-alpha));
ydot = @(x) (x);

% Parameters of numerical integration
Nstep = 5000;
h = Td/Nstep;

intervals = 100;
tmax = intervals*Td;

% Initalization
x_slide = (1:1:Nstep)*0.0;
y_slide = (1:1:Nstep)*0.0;
x_vals = (1:1:Nstep*intervals)*0.0;
y_vals = (1:1:Nstep*intervals)*0.0;
t_vals = (1:1:Nstep*intervals)*0.0;

t = 0.0;
init_cond = 0;
for i = 1:1:Nstep
   x_slide(i) = init_cond;
   y_slide(i) = init_cond;
end

% RK4 Method
i = 1;
while (t<tmax)
    
    k1x = xdot(x_slide(Nstep), x_slide(1), y_slide(Nstep));
    k1y = ydot(x_slide(Nstep));
    
    k2x = xdot(x_slide(Nstep)+h/2*k1x, x_slide(1), y_slide(Nstep)+h/2*k1y);
    k2y = ydot(x_slide(Nstep)+h/2*k1x);
    
    k3x = xdot(x_slide(Nstep)+h/2*k2x, x_slide(1), y_slide(Nstep)+h/2*k2y);
    k3y = ydot(x_slide(Nstep)+h/2*k2x);
    
    k4x = xdot(x_slide(Nstep)+h*k3x, x_slide(1), y_slide(Nstep)+h*k3y);
    k4y = ydot(x_slide(Nstep)+h*k3x);

    % Updating variables and storing in vectors.
    t = t + h;
    x_slide(Nstep) = x_slide(Nstep) + h/6*(k1x+2*k2x+2*k3x+k4x);
    y_slide(Nstep) = y_slide(Nstep) + h/6*(k1y+2*k2y+2*k3y+k4y);
    x_vals(i) = x_slide(Nstep);
    y_vals(i) = y_slide(Nstep);
    t_vals(i) = t;
    
    i = i + 1;
    
    % Sliding vector updates
    for j=1:1:(Nstep-1)
        x_slide(j) = x_slide(j+1);
        y_slide(j) = y_slide(j+1);
    end
    
end


% Plots
figure()
hold on
plot(x_vals, y_vals, 'LineWidth', 1.5);
plot(x_vals(1),y_vals(1),'mo','LineWidth', 1.5) % starting point
plot(x_vals(end),y_vals(end),'g*','LineWidth', 1.5) % ending point
legend(["RK4 Numerical Integration", "Starting Point", "Ending Point"], "location", "best")
title("Phase Portrait")
xlabel("x")
ylabel("y")
hold off

figure()
plot(t_vals, x_vals)
title("V_{RF}(t) Behavior")
xlabel("t")
ylabel("x(t)")

