%***** RUN 1D ADVECTION DIFFUSION MODEL ***********************************

% clear workspace
clear; close all;

% set CFL numbers to tune time step
CC = [1/1,1/2,1/4]/8;

for cc = 1:3

% set model parameters
W     = 1e3;           % domain width [m]
D     = 1e3;           % domain depth [m]
Nz    = 300;           % grid size z-direction
Nx    = Nz*W/D;        % grid size x-direction
h     = D/Nz;          % grid spacing (h = dx = dz)

kT0 = 2;               % thermal conductivity [W/m/K]
alphaT0 = 1e-6;        % thermal expansivity  [1/K]
rho0 = 2700;           % density [kg/m3]
cP0 = 1100;            % heat capacity [J/kg/K]
Qr0 = 1e-6;            % heat productivity [W/m3]
u0 = 1e-6;             % advection x-speed [m/s]
w0 = 1e-6;             % advection z-speed [m/s]


T0    = 100;           % initial background temperature [C]
dT    = 1000;          % initial temperature peak amplitude [C]
sgm0  = 25;            % initial temperature peak half-width (std dev.) [m]

k0    = kT0/(rho0*cP0);% heat diffusivity [m2/s]

BC     = 'periodic';   % boundary condition option flag ('insulating', 'periodic')
ADVN   = 'UPW3';       % advection scheme ('UPW1', 'CFD2', 'UPW3')
TINT   = 'RK2';        % time integration scheme ('FE1', 'RK2') (Implicit:'BE1', 'CN2')
SCHEME = 'explicit';   % Explicit or implicit scheme ('explicit', 'implicit')

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