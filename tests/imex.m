% Implicit and Explicit
% rk2+cranknicolson vs rk4

close all;
clear all;
dt = 0.0001;
N = 10;

y = zeros(1,N);
initial = 0.01;
y(1) = initial;

tmax = 100000;
time = 0:dt:tmax*dt;

for t = 1:tmax
    ystar = rk2(y(t),@ynlpt,dt);
    y(t+1) = ( (1+dt/2)/(1-dt/2) )* ystar;
end

z = initial;
z = rk4dyn(z,@dyn,dt,tmax);

figure;
hold on;
plot(time(1:1000:end),y(1:1000:end), 'o','DisplayName','Crank Nicolson');
plot(time,z, 'DisplayName','rk4', 'LineWidth',1);
legend;
title('dy/dt = y - y^3');

%% Functions

function y = ynlpt(y)
    y = - y^3;
end

function u = rk2(u,dynfunc,dt)
    k1 = dynfunc(u);
    u1 = u+k1*dt;
    k2 = dynfunc(u1);
    u = u + dt*((k1+k2)/2);
end

function dydt = dyn(y)
   dydt = y - y^3;
end

function z = rk4dyn(z,dynFunc,h,tmax)
% dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)
for t = 1:tmax
    % k values
    k1 = dynFunc(z(t));
    k2 = dynFunc(z(t) + (0.5*h)*k1);
    k3 = dynFunc(z(t) + (0.5*h)*k2);
    k4 = dynFunc(z(t) + h*k3);
    % dynamics
    z(t+1) = z(t) + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
end
end