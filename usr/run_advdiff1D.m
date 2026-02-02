%*****  RUN 2D ADVECTION DIFFUSION MODEL  *********************************

% clear workspace
clear; close all;

% set model parameters
W     = 1e3;           % domain width [m]
D     = 1e3;           % domain depth [m]
Nz    = 100;           % grid size z-direction
Nx    = Nz*W/D;        % grid size x-direction
h     = D/Nz;           % grid spacing (h = dx = dz)

kT0   = 2;             % thermal conductivity [W/m/K]

rho0  = 2700;          % density [kg/m3]
cP0   = 1100;          % heat capacity [J/kg/K]
Qr0   = 0e-6;          % heat productivity [W/m3]
u0    = 1e-6;          % advection x-speed [m/s]
w0    = 1e-6;          % advection z-speed [m/s]

% darcy flow parameters
g       = 9.81;        % gravity [m/s^2]
alphaT0 = 1e-6;        % thermal expansivity  [1/K]
KD0     = 1e-7;        % Darcy mobility [m^2/(Pa*s)]
Ttop    = 0;           % top temperature [C]
Tbot    = 1000;        % bottom temperature [C]
Tref    = Ttop;        % reference temperature for density law

T0    = 100;           % initial background temperature [C]

dT    = 1000;          % initial temperature peak amplitude [C]
sgm0  = 25;            % initial temperature peak half-width (std dev.) [m]

k0    = kT0/(rho0*cP0);% heat diffusivity[m2/s]

BC    = 'X_PerZ_Iso';    % boundary condition option flag ('insulating', 'periodic', 'X_PerZ_Iso)
ADVN  = 'UPW3';        % advection scheme ('UPW1', 'CFD2', 'UPW3')
TINT  = 'RK2';         % time integration scheme (Explicit:'FE1', 'RK2') (Implicit:'BE1', 'CN2')
SCHEME= 'explicit';    % Implicit or explicit scheme ('explicit', 'implicit')

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

