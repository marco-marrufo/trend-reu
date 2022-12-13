%%
% Author: Marco Marrufo
% This is a naive numerical simulation used to find the nth critical gain
% value that results in the nth nontrivial pair of fixed points for
% x/beta = cos^2[x+phi].
clear, clc

syms x

N_iter = 1000;
phi = -pi/2;
pair_numb = 1;

beta_init = 3;
beta = beta_init;
kbeta = 1.005;

start_lobe = pi + (pair_numb-1)*pi;
end_lobe = 3*pi/2 + (pair_numb-1)*pi;

lobe_range = [start_lobe end_lobe];

x1 = 0:0.01:6*pi;
f1 = cos(x1+phi).^2;

figure
hold on
plot(x1, f1)
title("Critical Gain Simulation")
xlabel("x")
ylabel("f(x)")
ylim([0, 1])
i = 0;
while i<N_iter
    disp("Iteration:" + i)
    eqn = x/beta == cos(x+phi)^2;
    solx = vpasolve(eqn,x,lobe_range);
    f2 = 1/beta*x1;
    plot(x1, f2, '--')
    if ~isempty(solx)
        disp("Solution Found! Beta (critical) = " + beta)
        break
    else
        i = i+1;
        beta = beta * kbeta;
    end
end