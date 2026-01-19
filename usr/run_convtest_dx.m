%***** RUN 1D ADVECTION DIFFUSION MODEL ***********************************

% clear workspace
clear all; close all; %clc;

NN = [200,400,800];

for nn = 1:3

% set model parameters
W     = 1000;          % domain width [m]
N     = NN(nn);        % grid size
dx    = W/N;           % grid spacing

T0    = 100;           % initial background temperature [C]
dT    = 1000;          % initial temperature peak amplitude [C]
sgm0  = 25;            % initial temperature peak half-width (std dev.) [m]

k0    = 0e-6;          % heat diffusivity [m2/s]
u0    = 1e-6;          % advection speed [m/s]

BC    = 'periodic';    % boundary condition option flag ('insulating', 'periodic')
ADVN  = 'UPW3';        % advection scheme ('UPW1', 'CFD2', 'UPW3')
TINT  = 'RK2';         % time integration scheme ('FE1', 'RK2')

yr    = 3600*24*365;   % seconds per year [s]
tend  = W/max(u0,k0);  % stopping time [s]
CFL   = N/10000;        % time step limiter
nop   = 5000;          % make output figure every 'nop' time steps

%*****  RUN MODEL
run('../src/main_solution.m');

E(nn)  = Err;
DX(nn) = dx;

end

figure(); 
loglog(DX,E                   ,'ro','LineWidth',2.0,'MarkerSize',8); axis tight; box on; hold on
loglog(DX,E(1).*[1,1/2,1/4].^1,'k-','LineWidth',1.5)
loglog(DX,E(1).*[1,1/2,1/4].^2,'k-','LineWidth',1.0)
loglog(DX,E(1).*[1,1/2,1/4].^3,'k-','LineWidth',0.5)
legend('num. error','linear','quadratic','cubic','FontSize',15,'box','off','location','southeast')
xlabel('Step size [m]','FontSize',18)
ylabel('Numerical error','FontSize',18)
title('Numerical Convergence in Space','FontSize',20)