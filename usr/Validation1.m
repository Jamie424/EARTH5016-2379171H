%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

% load model setup from image, reduce to target grid size
n_units = 10;                    % number of rock units contained in image
units   = imread('units.tiff');  % read in cross section with rock units
reduce  = 2;                     % reduce model resolution by factor
units   = imresize(units,1/reduce,'method','nearest');  % resize units map to model size

% set model dimensions
W       = 16e3;        % domain width (corresponds to width of geological cross section) [m]
[Nz,Nx] = size(units); % grid size according to units map
D       = W*Nz/Nx;     % domain depth (according to grid aspect ratio)
h       = W/Nx;        % grid spacing

% material properties for each rock unit (update based on your calibration)

matprop = [
% unit  kT    rho    cP     Qr    KD
   1   0.5   1000   1000   0e-6   1e-9   % air/water
   2   1.37	 2078	1625   0e-6   9e-8   % Sa
   3   1.43	 2153	1398   0e-6   4e-10  % Si
   4   1.93	 2106	1802   0e-6   1e-7   % Gr
   5	2.7	 2500	 820   4.5e-6 1e-11   % He1
   6	2.7	 3000	1000   0e-6   1e-11   % Bg
   7	2.7	 2600	830    6e-6   1e-11   % He2
   8    2.5	 2300	1000   0e-6   3.58e-7 % Fz
   9    1.7  2037	1451   0e-6   1e-9   % Ms
  10    1.9  2133	1209   0e-6   1e-10]; % Cm
air = units==1;
fault = units==8;


% get coefficient fields based on spatial distribution of rock units from image
kT0  = reshape(matprop(units,2),Nz,Nx);  % thermal conductivity [W/m/K]
rho0 = reshape(matprop(units,3),Nz,Nx);  % density [kg/m3]
cP0  = reshape(matprop(units,4),Nz,Nx);  % heat capacity [J/kg/K]
Qr0  = reshape(matprop(units,5),Nz,Nx);  % heat productivity [W/m3]
KD0  = reshape(matprop(units,6),Nz,Nx);  % segregation mobility [m2/Pas]
KD0  = imgaussfilt(KD0,1);               % apply some smoothing for numerical stability


u0    = 1e-6;          % advection x-speed [m/s]
w0    = 1e-6;          % advection z-speed [m/s]

% darcy flow parameters
g       = 9.81;        % gravity [m/s^2]
alphaT0 = 1e-5;        % thermal expansivity  [1/K]
rhoW0   = 1000;        % density of water [kg/m3]
Ttop    = 15;          % top temperature [C]
Tref    = Ttop;        % reference temperature for density law

BC    = 'insulating';  % boundary condition option flag ('insulating', 'periodic')
ADVN  = 'UPW3';        % advection scheme ('UPW1', 'CFD2', 'UPW3')
TINT  = 'RK2';         % time integration scheme (Explicit:'FE1', 'RK2')
MODE  = 'SIM';         % Verfification or simulation ('VERIFY','SIM')

yr    = 3600*24*365;   % seconds per year [s]
tend  = 1e7*yr;        % stopping time [s]
CFL   = 1/2;           % time step limiter
nop   = 20;            % make output figure every 'nop' time steps
tolP  = 1e-6;          % Pressure tolerance [Pa] 
alpha = 0.80;          % 
beta  = 0.80;
geo   = 0.03;          % geotherm [K/m]
qbot  = 0.1;          % vertical upward heat flux at boundary [W/m^2] 


% Diffusion and Basal Heat flux only (no radiogenic heating/Darcy
% mobility)   *****uncomment and run*****
% Qr0(:) = 0;    % No radiogenic heating
% KD0(:)   = 0;  % No Darcy flow
% alphaT0  = 0;  % No thermal expansivity
% u0 = 0; w0 = 0;



% drillhole data D
zdrill = [51; 150; 249; 351; 450; 550; 649; 749; 850; 950; 1049; 1150; 1253; 1352];
Tdrill = [ 15.1; 15.9; 21.5; 21.0; 25.0; 28.3; 31.6; 35.7; 40.4; 43.2; 45.6; 50.8; 55.4; 56.4];
x_dh = 5000;


% proposed drill site
x_pds = 11200;

% Turn on Darcy flow after x years
tDarcyOn = 100000 * yr; % [s]


runID = 'Conduction';
%*****  RUN MODEL
run('../src/main.m');