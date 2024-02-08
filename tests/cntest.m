% Crank Nicolson test with dy/dt = -y
% exact solution: y = exp(-t)

close all;
clear all;
dt = 0.01;
N = 10;

y = zeros(1,N);
y(1,1) = 1;

tmax = 1000;
time = 0:dt:tmax*dt;

for t = 1:tmax
    y(t+1) = ( (1-dt/2)/(1+dt/2) )* y(t);
end

figure;
plot(time(1:10:end),y(1:10:end), 'o','DisplayName','Crank Nicolson');
hold on;
plot(time,exp(-time), 'DisplayName','Exact', 'LineWidth',1);
legend;
title('dy/dt = -y');