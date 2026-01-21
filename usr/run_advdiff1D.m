%*****  RUN 1D ADVECTION DIFFUSION MODEL  *********************************

% clear workspace
clear; close all;

% set model parameters
W     = 1000;          % domain width [m]
N     = 400;           % grid size
dx    = W/N;           % grid spacing

T0    = 100;           % initial background temperature [C]
dT    = 1000;          % initial temperature peak amplitude [C]
sgm0  = 25;            % initial temperature peak half-width (std dev.) [m]

k0    = 1e-6;          % heat diffusivity [m2/s]
u0    = 1e-6;          % advection speed [m/s]

BC    = 'periodic';    % boundary condition option flag ('insulating', 'periodic')
ADVN  = 'UPW3';        % advection scheme ('UPW1', 'CFD2', 'UPW3')
TINT  = 'RK2';         % time integration scheme (Explicit:'FE1', 'RK2') (Implicit:'BE1', 'CE2')
SCHEME= 'implicit';   % Implicit or explicit scheme ('explicit', 'implicit')

yr    = 3600*24*365;   % seconds per year [s]
tend  = W/max(u0,k0);  % stopping time [s]
CFL   = 1/8;           % time step limiter
nop   = 100;           % make output figure every 'nop' time steps


%*****  RUN MODEL
run('../src/main.m');


%%%%%%%%%%%%%%%%%%%%%%

% Qualitive test for performance of ADVN schemes

% k0 = 0 no diffusion, u0 = 1e-6 advection speed, 
% periodic BC, RK2 time integration

% UPW1 suffers from heavy diffusion at longer time frames as there is
% smoothing and decay compared to analytical solution.

% centred differencing on the other hand is plagued by numerical dispersion
% as there are oscillations on the left, as the velocity is directional to
% the right and hence not symmetric

% UPW3 3rd order upwind scheme reduces the error rate and is noticeably
% closer to analytical solution. Higher accuracy as rate of change is
% evaluated and updated three times, decreasing the discretisation error 
% each time 

%%%%%%%%%%%%%%%%%%%%%%%%

% Errors to fix - the analytical solution dotted plot is not moving to the
% right with advection. It does however diffuse when adding a diffusion
% term. 

