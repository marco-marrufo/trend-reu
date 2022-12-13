%% RK4 for DDE Examples
% Code taken from Professor Chembo's lecture slides.
close, clc

% System parameters
gamma = 2.0;
phi = pi/2;
Td = 2.0;

% Parameters of numerical integration
Nstep = 2000;
h = Td/Nstep;

intervals = 100;
tmax = intervals*Td;


% Initalization
x = (1:1:Nstep)*0.0;
x_vals = (1:1:Nstep*intervals)*0.0;
t_vals = (1:1:Nstep*intervals)*0.0;

t = 0.0;
for i = 1:1:Nstep
   x(i) = 0.01; 
end

figure()
hold on
% RK4 Method
i = 1;
while (t<tmax)
    
    k1 = -x(Nstep) + gamma * sin(x(1) + phi);
    k2 = -(x(Nstep)+0.5*h*k1) + gamma * sin(x(1) + phi);
    k3 = -(x(Nstep)+0.5*h*k2) + gamma * sin(x(1) + phi);
    k4 = -(x(Nstep)+h*k3) + gamma * sin(x(1) + phi);
    
    k = (k1+2*k2+2*k3+k4)/6;
    x(Nstep) = x(Nstep) + k*h;
    t = t+h;
    x_vals(i) = x(Nstep);
    t_vals(i) = t;
    
    i = i + 1;
    
    % Sliding vector updates
    for j=1:1:(Nstep-1)
        x(j) = x(j+1);
    end
    
end

plot(t_vals, x_vals)
