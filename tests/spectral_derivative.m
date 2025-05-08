% spectral derivative
close all;
clear all;

L = 20;
N = 100;

x = linspace(-L/2,L/2,N+1);
x = x(1:N);
kx = (2*pi/L)*(-N/2:N/2-1);
kx = fftshift(kx);

%% Periodic function
f = cos((2*pi/L)*x);

fdash = 1i*kx.*fft(f);
fdash = real(ifft(fdash));

figure; 
subplot(2,1,1);
plot(x,fdash,'-^','DisplayName','df/dx (spectral)','LineWidth',1.5);
hold on; plot(x,-(2*pi/L)*sin((2*pi/L)*x), ...
    'DisplayName','-(2*pi/L)*sin((2*pi/L)*x)',LineWidth=1.5);
legend;
set(gca,'TickLabelInterpreter','tex','FontSize',15);
xlabel('$x$', 'Interpreter','latex','FontSize',30);
title('$df/dx \, for \, f(x) = \cos(2\pi x/L) \, (periodic)$', ...
    'Interpreter','latex');

%% Aperiodic function
f2 = cos((3*pi/L)*x);
f2dash = 1i*kx.*fft(f2);
f2dash = real(ifft(f2dash));

subplot(2,1,2);
plot(x,f2dash,'-^','DisplayName','df/dx (spectral)','LineWidth',1.5);
hold on; plot(x,-(3*pi/L)*sin((3*pi/L)*x), ...
    'DisplayName','-(2*pi/L)*sin((2*pi/L)*x)',LineWidth=1.5);
legend;
set(gca,'TickLabelInterpreter','tex','FontSize',15);
xlabel('$x$', 'Interpreter','latex','FontSize',30);
title('$df/dx \, for \, f(x) = \cos(3\pi x/L) \, (aperiodic)$', ...
    'Interpreter','latex');