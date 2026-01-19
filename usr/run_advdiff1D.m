%*****  RUN 1D ADVECTION DIFFUSION MODEL  *********************************

% clear workspace
clear; close all;

% set model parameters
W     = 1000;          % domain width [m]
N     = 200;           % grid size
dx    = W/N;           % grid spacing

T0    = 100;           % initial background temperature [C]
dT    = 1000;          % initial temperature peak amplitude [C]
sgm0  = 25;            % initial temperature peak half-width (std dev.) [m]

k0    = 0e-6;          % heat diffusivity [m2/s]
u0    = 1e-6;          % advection speed [m/s]

BC    = 'periodic';    % boundary condition option flag ('insulating', 'periodic')
ADVN  = 'UPW3';        % advection scheme ('UPW1', 'CFD2', 'UPW3')
TINT  = 'FE1';         % time integration scheme ('FE1', 'RK2')

yr    = 3600*24*365;   % seconds per year [s]
tend  = W/max(u0,k0)/4;  % stopping time [s]
CFL   = 1/4;           % time step limiter
nop   = 100;           % make output figure every 'nop' time steps


%*****  RUN MODEL
run('../src/main_solution.m');