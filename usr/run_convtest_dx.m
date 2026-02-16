%***** RUN 1D ADVECTION DIFFUSION MODEL ***********************************

% clear workspace
clear all; close all; %clc;

NN = [200,400,800];

for nn = 1:3

% set model parameters
W     = 1e3;           % domain width [m]
D     = 1e3;           % domain depth [m]
Nz    = NN(nn);        % grid size z-direction
Nx    = Nz*W/D;        % grid size x-direction
h     = D/Nz;          % grid spacing (h = dx = dz)

kT0   = 2;             % thermal conductivity [W/m/K]
alphaT0 = 1e-6;        % thermal expansivity  [1/K]
rho0  = 2700;          % density [kg/m3]
cP0   = 1100;          % heat capacity [J/kg/K]
Qr0   = 0e-6;          % heat productivity [W/m3]
u0    = 1e-6;          % advection x-speed [m/s]
w0    = 1e-6;          % advection z-speed [m/s]

KD0     = 1e-7;        % Darcy mobility [m^2/(Pa*s)]

T0    = 100;           % initial background temperature [C]
dT    = 1000;          % initial temperature peak amplitude [C]
sgm0  = 25;            % initial temperature peak half-width (std dev.) [m]

k0    = kT0/(rho0*cP0);% heat diffusivity[m2/s]

BC    = 'periodic';    % boundary condition option flag ('insulating', 'periodic')
ADVN  = 'UPW3';        % advection scheme ('UPW1', 'CFD2', 'UPW3')
TINT  = 'RK2';         % time integration scheme ('FE1', 'RK2')
SCHEME= 'explicit';    % Explicit or implicit scheme ('explicit', 'implicit')
MODE   = 'VERIFY';     % Verfification ('VERIFY','SIM')

yr    = 3600*24*365;   % seconds per year [s]
tend  = 0.1 * W/max(u0,k0);  % stopping time [s]
CFL   = 1/10;         % time step limiter
nop   = 5000;          % make output figure every 'nop' time steps

%*****  RUN MODEL
run('../src/main.m');

E(nn)  = Err;
DX(nn) = h;

end

DXref = DX(1) ./ [1,2,4]; 
figure(); 
loglog(DX,E                   ,'ro','LineWidth',2.0,'MarkerSize',8); axis tight; box on; hold on  
loglog(DXref, E(1)*(DXref/DX(1)).^1, 'k-', 'LineWidth',1.5)
loglog(DXref, E(1)*(DXref/DX(1)).^2, 'k-', 'LineWidth',1.0)
loglog(DXref, E(1)*(DXref/DX(1)).^3, 'k-', 'LineWidth',0.5)
% loglog(DX,E(1).*[200,400,800].^1,'k-','LineWidth',1.5)
% loglog(DX,E(1).*[200,400,800].^2,'k-','LineWidth',1.0)
% loglog(DX,E(1).*[200,400,800].^3,'k-','LineWidth',0.5)
legend('num. error','linear','quadratic','cubic','FontSize',15,'box','off','location','southeast')
xlabel('Step size [m]','FontSize',18)
ylabel('Numerical error','FontSize',18)
title('Numerical Convergence in Space','FontSize',20)