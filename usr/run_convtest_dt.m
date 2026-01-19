%***** RUN 1D ADVECTION DIFFUSION MODEL ***********************************

% clear workspace
clear; close all;

% set CFL numbers to tune time step
CC = [1/1,1/2,1/4]/8;

for cc = 1:3

% set model parameters
W     = 1000;          % domain width [m]
N     = 4000;           % grid size
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
CFL   = CC(cc);        % time step limiter
nop   = 1000;          % make output figure every 'nop' time steps

%*****  RUN MODEL
run('../src/main.m');

E(cc)  = Err;
DT(cc) = dt;

end

figure(); 
loglog(DT,E                   ,'ro','LineWidth',2.0,'MarkerSize',8); axis tight; box on; hold on
loglog(DT,E(1).*[1,1/2,1/4].^1,'k-','LineWidth',1.5)
loglog(DT,E(1).*[1,1/2,1/4].^2,'k-','LineWidth',1.0)
loglog(DT,E(1).*[1,1/2,1/4].^3,'k-','LineWidth',0.5)
legend('num. error','linear','quadratic','cubic','FontSize',15,'box','off','location','southeast')
xlabel('Step size [s]','FontSize',18)
ylabel('Numerical error','FontSize',18)
title('Numerical Convergence in Time','FontSize',20)