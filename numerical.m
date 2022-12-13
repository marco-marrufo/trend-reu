%% Numerical Methods Practice
clear, clc

dydx = @(x, y) 3*x^2*y;

% Step size
h = 0.01;

% Initial conditions
x0 = 1;
y0 = 2;


% X max
xmax = 2;
iters = int64((xmax - x0)/h)+1;

%% Euler's Method

disp("Running Euler's Method")

x = x0;
y = y0;


euler_xvals = zeros(1, iters);
euler_yvals = zeros(1, iters);

i = 1;

euler_xvals(i) = x0;
euler_yvals(i) = y0;

while x < xmax
    i = i + 1;
    
    x_prev = x;
    y_prev = y;
    
    x = x_prev + h;
    y = y_prev + h*dydx(x_prev, y_prev);
    
    euler_xvals(i) = x;
    euler_yvals(i) = y;
    
end


%% RK4

disp("Running RK4")

x = x0;
y = y0;

rk4_xvals = zeros(1, iters);
rk4_yvals = zeros(1, iters);

i = 1;

rk4_xvals(i) = x0;
rk4_yvals(i) = y0;

while x < xmax
    i = i + 1;
    
    x_prev = x;
    y_prev = y;
    
    % RK4 Calculations
    k1 = dydx(x_prev, y_prev);
    k2 = dydx(x_prev + h/2, y_prev + h/2*k1);
    k3 = dydx(x_prev + h/2, y_prev + h/2*k2);
    k4 = dydx(x_prev + h, y_prev + h*k3);
    
    x = x_prev + h;
    y = y_prev + h/6*(k1+2*k2+2*k3+k4);
    
    rk4_xvals(i) = x;
    rk4_yvals(i) = y;
end

%%
figure
plot(euler_xvals, euler_yvals, 'k.');
hold on
plot(rk4_xvals, rk4_yvals);
title("Numerical Methods for $$\frac{dy}{dx} = 3x^2y$$", 'interpreter','latex')
legend(["Euler's", "RK4"], "Location", "northwest")
xlim([x0-h, xmax+h])
xlabel("x")
ylabel("y")