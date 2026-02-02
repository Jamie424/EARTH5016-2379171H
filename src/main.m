%*****  2D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT  *******************


%*****  Initialise Model Setup

% create coordinate vectors
xc = h/2:h:W-h/2;          % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;          % z-coordinate vector for cell centre positions [m]
xf = 0:h:W ;               % x-coordinate vector for cell face positions [m]
zf = 0:h:D ;               % z-coordinate vector for cell face positions [m]
[Xc,Zc] = meshgrid(xc,zc); % create 2D coordinate arrays



% set up ghosted index lists for boundary conditions
switch BC
    case 'periodic'
        % 3-point stencil            |-- i-1 --|-- i --|-- i+1 --|
        % 5-point stencil  |-- i-2 --|-- i-1 --|-- i --|-- i+1 --|-- i+2 --|
        ix3 = [      Nx, 1:Nx, 1   ];  iz3 = [      Nz, 1:Nz, 1   ];
        ix5 = [Nx-1, Nx, 1:Nx, 1, 2];  iz5 = [Nz-1, Nz, 1:Nz, 1, 2];
        
    case 'insulating'
        % example non-periodic indexing for N=4 
        % 3-point stencil            |-- i-1 --|-- i --|-- i+1 --|
        % 5-point stencil  |-- i-2 --|-- i-1 --|-- i --|-- i+1 --|-- i+2 --|
        ix3 = [   1, 1:Nx, Nx    ];    iz3 = [   1, 1:Nz, Nz    ];
        ix5 = [1, 1, 1:Nx, Nx, Nx];    iz5 = [1, 1, 1:Nz, Nz, Nz]; 

    case 'X_PerZ_Iso'
        % Isothermal at top/bottom, periodic horizontally
        ix3 = [      Nx, 1:Nx, 1   ];    iz3 = [   1, 1:Nz, Nz    ];
        ix5 = [Nx-1, Nx, 1:Nx, 1, 2];    iz5 = [1, 1, 1:Nz, Nz, Nz]; 
end



% set initial coefficient fields
kT     = kT0     .* ones(Nz,Nx);         % matrix of initial thermal conductivity [W/m/K]
alphaT = alphaT0 .* ones(Nz,Nx);         % matrix of initial thermal expansivity  [W/m/K]
cP     = cP0     .* ones(Nz,Nx);         % heat capacity [J/kg/K]
rho    = rho0    .* ones(Nz,Nx);         % density [kg/m3]
Qr     = Qr0     .* ones(Nz,Nx);         % heat production [W/m3]

% pressure field
p = zeros(Nz,Nx); 

% KD Darcy mobility field
KD = KD0 .* ones(Nz,Nx);


% Initial Darcy flux
u = zeros(Nz,Nx+1);          % x-face flux  
w = zeros(Nz+1,Nx);          % z-face flux



% set time step size (Courant–Friedrichs–Lewy condition)
dt_adv = (h/2)   / (max(abs(u(:))) + max(abs(w(:))) + eps); 
dt_dff = (h/2)^2 / max(kT(:)./rho(:)./cP(:) + eps);                                              
dt     = CFL * min(dt_adv,dt_dff);               % time step [s]

% Initialise time count variables
t = 0;  % initial time [s]
k = 0;  % initial time step count

% set initial temperature field (linear gradient from Tbot to Ttop)
T = Ttop + (zc(:)/D) * (Tbot - Ttop) .* ones(1,Nx);
T(1,:)   = Ttop;
T(end,:) = Tbot;
Tin = T;
Ta = Tin;    % dummy for now


% initialise output figure with initial condition
figure(1); clf
makefig(xc,zc,T,Tin,Ta,0);

%*****  Solve Model Equations
while t <= tend

    % increment time and step count
    t = t + dt;
    k = k + 1;

    % Temperature dependant density
    rho = rho0 .* (1 - alphaT.*(T - Tref));

    [u,w] = darcy_flux(p, rho, KD, h, g, ix3, iz3)

    % Closed top/bottom boundaries (no flow through vertical boundaries)
    w(1,:)   = 0;
    w(end,:) = 0;
    
   
    switch TINT % select explicit time integration scheme
        case 'FE1'  % 1st-order Forward Euler time integration scheme
            % get rate of change
            dTdt = diffusion(T,kT,h,ix3,iz3)./(rho.*cP)  ...
                 + advection(T,u,w,h,ix5,iz5,ADVN) ...
                 + Qr./(rho.*cP);
                  

        case 'RK2'  % 2nd-order Runge-Kutta time integration scheme
            dTdt_half = diffusion(T                ,kT,h,ix3,iz3)./(rho.*cP) ...
                      + advection(T               ,u,w,h,ix5,iz5,ADVN) ...
                      + Qr./(rho.*cP);
            dTdt      = diffusion(T+dTdt_half*dt/2,kT, h,ix3,iz3)./(rho.*cP) ...
                      + advection(T+dTdt_half*dt/2,u,w,h,ix5,iz5,ADVN) ...
                      + Qr./(rho.*cP);
                      
    end

    % update temperature
    T = T + dTdt * dt;  

    % update boundary conditions for temperature
    T(1,:)   = Ttop;
    T(end,:) = Tbot;

    % get analytical solution at time t
    Ta = analytical(T0,dT,sgm0,k0,u0,w0,Xc,Zc,D,W,t);

    % plot model progress
    if ~mod(k,nop)
        makefig(xc,zc,T,Tin,Ta,t/yr);
        pause(0.1);
    end


end

%*****  calculate and display numerical error norm
%Err = rms(T - Ta, 'all') / rms(Ta, 'all');

disp(' ');
disp(['Advection scheme: ',ADVN]);
disp(['Time integration scheme: ',TINT]);
disp(['Numerical error norm = ',num2str(Err)]);
disp(' ');



%*****  Utility Functions  ************************************************

%*****  Function to make output figure

function makefig(x,z,T,Tin,Ta,t)

subplot(2,1,1);
imagesc(x,z,T); axis equal tight; colorbar
ylabel('z [m]','FontSize',15)
title(['Temperature [C]; time = ',num2str(t),' [yr]'],'FontSize',17)
subplot(2,1,2)
imagesc(x,z,T - Ta); axis equal tight; colorbar
xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('Normalized Num. Error [C]','FontSize',17)

drawnow;

end

%**************************************************************************
%*****  Function to calculate diffusion rate

function dfdt = diffusion(f,k,h,ix,iz)

% input arguments
% f:    diffusing scalar field
% k:    diffusion coefficient
% h:   grid spacing
% ix:  x ghosted index list
% iz:  z ghosted index list

% output variables
% dfdt: diffusion rate of scalar field f

% calculate diffusive flux coefficient at cell faces
kfz = (k(iz(1:end-1),:)+k(iz(2:end),:))/2;
kfx = (k(:,ix(1:end-1))+k(:,ix(2:end)))/2;

% calculate diffusive flux of scalar field f
qz = - kfz .* diff(f(iz,:),1,1)/h;
qx = - kfx .* diff(f(:,ix),1,2)/h;

% calculate diffusion flux balance for rate of change
dfdt = - diff(qz,1,1)/h ...
       - diff(qx,1,2)/h;

% Note: start of function defines the output dfdt
end

%**************************************************************************
%*****  Function to calculate advection rate

function dfdt = advection(f,u,w,h,ix,iz,ADVN)

% input arguments
% f:    advected scalar field
% u:    advection velocity x direction 
% w:    advection velocity z direction
% h:    grid spacing
% ix:  ghosted index list x direction
% iz:  ghosted index list z direction
% ADVN: advection scheme ('UPW1', 'CFD2', 'UPW3')

% output variables
% dfdt: advection rate of scalar field


%****** horizontal advection (x direction)

% split the velocities into positive and negative
u_pos = max(0,u);    % positive velocity (to the right)
u_neg = min(0,u);    % negative velocity (to the left)

% get values on stencil nodes
f_imm  = f(:,ix(1:end-4));  % i-2
f_im   = f(:,ix(2:end-3));  % i-1   
f_ic   = f(:,ix(3:end-2));  % i
f_ip   = f(:,ix(4:end-1));  % i+1
f_ipp  = f(:,ix(5:end  ));  % i+2

% get interpolated field values on i+1/2, i-1/2 cell faces
switch ADVN
    case 'UPW1'   % 1st-order upwind scheme
        % positive velocity          -> boundary inherited from left
        f_ip_pos = f_ic;     % i+1/2 value of right face look at center i 
        f_im_pos = f_im;     % i-1/2 value of left face look at left node i-1
        
        % negative velocity          <- boundary inherited from right
        f_ip_neg = f_ip;     % i+1/2 value of right face look at right node i+1 
        f_im_neg = f_ic;     % i-1/2 value of left face look at center node i

    case 'CFD2'  % 2nd-order centred finite-difference scheme
        % positive velocity
        f_ip_pos = (f_ic+f_ip)./2;     % i+1/2
        f_im_pos = (f_ic+f_im)./2;     % i-1/2

        % negative velocity
        f_ip_neg = f_ip_pos;           % i+1/2
        f_im_neg = f_im_pos;           % i-1/2

    case 'UPW3'  % 3rd-order upwind scheme
        % positive velocity
        f_ip_pos = (2*f_ip + 5*f_ic - f_im )./6;     % i+1/2                Note: f_ip_pos - f_im_pos later to get the 3rd order upwind scheme in slides
        f_im_pos = (2*f_ic + 5*f_im - f_imm)./6;     % i-1/2     

        % negative velocity
        f_ip_neg = (2*f_ic + 5*f_ip - f_ipp)./6;     % i+1/2
        f_im_neg = (2*f_im + 5*f_ic - f_ip )./6;     % i-1/2
 end

% calculate advection fluxes on i+1/2, i-1/2 cell faces

% positive velocity
q_ip_pos = u_pos(:,2:end).*f_ip_pos;        % flux on right face i+1/2 
q_im_pos = u_pos(:,1:end-1).*f_im_pos;        % flux on left face  i-1/2

% negative velocity
q_ip_neg = u_neg(:,2:end).*f_ip_neg;        % flux on right face i+1/2 
q_im_neg = u_neg(:,1:end-1).*f_im_neg;        % flux on left face  i-1/2

% advection flux balance for rate of change
div_qx_pos = (q_ip_pos - q_im_pos)/h;  % positive velocity                  
div_qx_neg = (q_ip_neg - q_im_neg)/h;  % negative velocity

div_qx     = div_qx_pos + div_qx_neg;     % combined flux x direction


%****** vertical advection (z direction)


% split the velocities into positive and negative
w_pos = max(0,w);    % positive velocity ('up' direction)
w_neg = min(0,w);    % negative velocity ('down' direction)

f_jpp = f(iz(5:end  ),:);   % j + 2
f_jp  = f(iz(4:end-1),:);   % j + 1
f_jc  = f(iz(3:end-2),:);   % j
f_jm  = f(iz(2:end-3),:);   % j - 1
f_jmm = f(iz(1:end-4),:);   % j - 2

switch ADVN
    case 'UPW1'   % 1st-order upwind scheme
        % positive velocity          boundary inherited from below
        f_jp_pos = f_jc;     % j+1/2 face above center i
        f_jm_pos = f_jm;     % j-1/2 value of left face look at left node i-1
        
        % negative velocity          boundary inherited from above
        f_jp_neg = f_jp;     % j+1/2
        f_jm_neg = f_jc;     % j-1/2 

    case 'CFD2'  % 2nd-order centred finite-difference scheme
        % positive velocity
        f_jp_pos = (f_jc+f_jp)./2;     % j+1/2
        f_jm_pos = (f_jc+f_jm)./2;     % j-1/2

        % negative velocity
        f_jp_neg = f_jp_pos;           % j+1/2
        f_jm_neg = f_jm_pos;           % j-1/2

    case 'UPW3'  % 3rd-order upwind scheme
        % positive velocity
        f_jp_pos = (2*f_jp + 5*f_jc - f_jm )./6;     % j+1/2                Note: f_ip_pos - f_im_pos later to get the 3rd order upwind scheme in slides
        f_jm_pos = (2*f_jc + 5*f_jm - f_jmm)./6;     % j-1/2     

        % negative velocity
        f_jp_neg = (2*f_jc + 5*f_jp - f_jpp)./6;     % j+1/2
        f_jm_neg = (2*f_jm + 5*f_jc - f_jp )./6;     % j-1/2

end
% positive velocity
q_jp_pos = w_pos(2:end,:) .* f_jp_pos;      % flux on top    face j+1/2
q_jm_pos = w_pos(1:end-1,:) .* f_jm_pos;      % flux on bottom face j-1/2

% negative velocity
q_jp_neg = w_neg(2:end,:) .* f_jp_neg;      % flux on top    face j+1/2
q_jm_neg = w_neg(1:end-1,:) .* f_jm_neg;      % flux on bottom face j-1/2

% advection flux balance for rate of change (z-direction)
div_qz_pos = (q_jp_pos - q_jm_pos)/h;   % positive velocity
div_qz_neg = (q_jp_neg - q_jm_neg)/h;   % negative velocity

div_qz     = div_qz_pos + div_qz_neg;   % combined flux x direction

div_q = div_qx + div_qz;                % x flux + z flux
dfdt  = - div_q;            % advection rate
end


% Function for Darcy flux (copilot help)
function [u,w] = darcy_flux(p, rho, KD, h, g, ix3, iz3)

% Computes Darcy flux components vD = [u,w]

% Darcy law used:
%   u = -KD * dp/dx
%   w = -KD * (dp/dz - rho*g)

% INPUTS:
%   p   : (Nz x Nx) pressure [Pa] (cell centres)
%   rho : (Nz x Nx) density  [kg/m^3] (cell centres)
%   KD  : (Nz x Nx) mobility [m^2/(Pa*s)] (cell centres)
%   h   : grid spacing [m]
%   g   : gravity [m/s^2] (positive, acts downward)
%   ix3 : ghosted x-index list length Nx+2 (periodic)
%   iz3 : ghosted z-index list length Nz+2 (choose for BC)

% OUTPUTS:
%   u : (Nz x (Nx+1)) x-face Darcy flux [m/s]
%   w : ((Nz+1) x Nx) z-face Darcy flux [m/s]

    % x-faces pressures:
    pL = p(:, ix3(1:end-1));     % Nz x (Nx+1)
    pR = p(:, ix3(2:end  ));     % Nz x (Nx+1)

    dpdx = (pR - pL) / h;        

    % face mobility: average KD from adjacent cells
    KDx = 0.5 * (KD(:, ix3(1:end-1)) + KD(:, ix3(2:end)));

    % Darcy flux in x
    u = -KDx .* dpdx;


    % z-faces
    % bottom/top cell pressures around each z-face
    pB = p(iz3(1:end-1), :);     % (Nz+1) x Nx
    pT = p(iz3(2:end  ), :);     % (Nz+1) x Nx
    dpdz = (pT - pB) / h;        

    % face mobility and density: average from adjacent cells
    KDz  = 0.5 * (KD(iz3(1:end-1), :)  + KD(iz3(2:end), :));
    rhoz = 0.5 * (rho(iz3(1:end-1), :) + rho(iz3(2:end), :));

    % Darcy flux in z (includes buoyancy)
    w = -KDz .* (dpdz - rhoz * g);
end

function F = darcy_residual(u,w,h)
% u: Nz x (Nx+1)
% w: (Nz+1) x Nx
% F: Nz x Nx (divergence at cell centres)
    F = diff(w)/h + diff(u)/h;
end


% function below not used for Darcy flux version - no known analytical solution

%****function to calculate the analytical temperature at time t 
function Ta = analytical(f0,df,sgm0,k0,u0,w0,Xc,Zc,D,W,t)

sgmt = sqrt(sgm0^2 + 2*k0*t);

% sum each of the 9 combinations
Ta   = f0 + df*(sgm0.^2/sgmt.^2)*(exp(-((Xc-(W/2    )- u0*t).^2 + (Zc-(D/2    ) - w0*t).^2)./ (2*sgmt^2)) ...
                                + exp(-((Xc-(W/2 + W)- u0*t).^2 + (Zc-(D/2    ) - w0*t).^2)./ (2*sgmt^2)) ...
                                + exp(-((Xc-(W/2 - W)- u0*t).^2 + (Zc-(D/2    ) - w0*t).^2)./ (2*sgmt^2)) ...
                                + exp(-((Xc-(W/2    )- u0*t).^2 + (Zc-(D/2 + D) - w0*t).^2)./ (2*sgmt^2)) ...
                                + exp(-((Xc-(W/2    )- u0*t).^2 + (Zc-(D/2 - D) - w0*t).^2)./ (2*sgmt^2)) ...
                                + exp(-((Xc-(W/2 + W)- u0*t).^2 + (Zc-(D/2 + D) - w0*t).^2)./ (2*sgmt^2)) ...
                                + exp(-((Xc-(W/2 + W)- u0*t).^2 + (Zc-(D/2 - D) - w0*t).^2)./ (2*sgmt^2)) ...
                                + exp(-((Xc-(W/2 - W)- u0*t).^2 + (Zc-(D/2 + D) - w0*t).^2)./ (2*sgmt^2)) ...
                                + exp(-((Xc-(W/2 - W)- u0*t).^2 + (Zc-(D/2 - D) - w0*t).^2)./ (2*sgmt^2)));

end

