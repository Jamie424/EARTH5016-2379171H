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
end

% set initial coefficient fields
kT     = kT0     .* ones(Nz,Nx);         % matrix of initial thermal conductivity [W/m/K]
alphaT = alphaT0 .* ones(Nz,Nx);         % matrix of initial thermal expansivity  [W/m/K]
cP     = cP0     .* ones(Nz,Nx);         % heat capacity [J/kg/K]
rho    = rho0    .* ones(Nz,Nx);         % density [kg/m3]
Qr     = Qr0     .* ones(Nz,Nx);         % heat production [W/m3]

% pressure field
p = zeros(Nz,Nx); 
resP = p*0;

switch MODE
    case 'VERIFY'
        KD = zeros(Nz,Nx);
        qbot = 0;
        air = false(Nz,Nx);

        % constant velocities
        w = w0 * ones(Nz+1, Nx);
        u = u0 * ones(Nz, Nx+1);


    case 'SIM'
        % KD Darcy mobility field
        KD = KD0 .* ones(Nz,Nx);
        % Initial Darcy flux
        u = zeros(Nz,Nx+1);          % x-face flux  
        w = zeros(Nz+1,Nx);          % z-face flux
end

% set initial temperature field
switch MODE
    case 'VERIFY'
        Ta = analytical(T0,dT,sgm0,k0,u0,w0,Xc,Zc,D,W,0);
        T=Ta;
    case 'SIM'
        dr = randn(Nz,Nx);
        dr = imgaussfilt(dr,1);
        T = Ttop + geo * (zc(:) * ones(1,Nx)) + dr*1;
        T(air) = Ttop;
end
Tin  = T;
resT = T*0;
dTdt = 0*T;

% initialise darcy velocity
switch MODE
    case 'SIM'
    rhoW = rhoW0 .* (1 - alphaT.*(T - Tref));
    Drho = rhoW - mean(rhoW,2);
    [u,w,resP] = darcy_flux(p, Drho, g, KD, h, ix3, iz3);
end

% set time step size (Courant–Friedrichs–Lewy condition)
dt_adv = (h/2)   / (max(abs(u(:))) + max(abs(w(:))) + eps); 
dt_dff = (h/2)^2 / max(kT(:)./rho(:)./cP(:) + eps);                                              
dt     = CFL * min(dt_adv,dt_dff);               % time step [s]

dt_init = dt;

% pseudo-transient iterative step
dtau = (h/2)^2 ./ max(KD(:) + eps);

% Initialise time count variables
t = 0;  % initial time [s]
k = 0;  % initial time step count


switch MODE
    case 'VERIFY'
        % initialise output figure with initial condition
        figure(1); clf
        makefigverify(xc,zc,T,Ta,t/yr);
    case 'SIM'
        figure(2); clf
        initialfig(xc,zc,T,KD,Qr,kT,t/yr);
        
end


%*****  Solve Model Equations
while t <= tend
    
    % Print time step information
    fprintf(1,'*****  step %d;  dt = %4.4e;  time = %4.4e [yr]\n\n',k,dt/yr,t/yr)

    % store previous T and dTdt
    To = T;
    k = k + 1; 
    
    switch MODE

        % constant velocity, No Darcy flow, pure advection diffusion
        case 'VERIFY'
                switch TINT % select explicit time integration scheme
                    case 'FE1'  % 1st-order Forward Euler time integration scheme
                        % get rate of change
                        dTdt = diffusion(T,kT,h,ix3,iz3,qbot)./(rho.*cP)  ...
                             + advection(T,u,w,h,ix5,iz5,ADVN) ...
                             + Qr./(rho.*cP);
            
            
                    case 'RK2'  % 2nd-order Runge-Kutta time integration scheme
                        dTdt_half = diffusion(T                ,kT,h,ix3,iz3,qbot)./(rho.*cP) ...
                                  + advection(T               ,u,w,h,ix5,iz5,ADVN) ...
                                  + Qr./(rho.*cP);
                        dTdt      = diffusion(T+dTdt_half*dt/2,kT, h,ix3,iz3,qbot)./(rho.*cP) ...
                                  + advection(T+dTdt_half*dt/2,u,w,h,ix5,iz5,ADVN) ...
                                  + Qr./(rho.*cP);
                end
                % update temperature
                T = T + dTdt * dt;  

                % *****  calculate and display numerical error norm
                Ta  = analytical(T0,dT,sgm0,k0,u0,w0,Xc,Zc,D,W,t);
                Err = rms(T - Ta, 'all') / rms(Ta, 'all');
                
                disp(' ');
                disp(['Advection scheme: ',ADVN]);
                disp(['Time integration scheme: ',TINT]);
                disp(['Numerical error norm = ',num2str(Err)]);
                disp(' ');
               

                if ~mod(k,nop)
                    makefigverify(xc,zc,T,Ta,t/yr)
                    pause(0.1);
                end
                         

        % Full Darcy/geothermal model
        case 'SIM'

            %%%%%%%%%%%%--------------
            % if time < x years then Darcy is off
            %%%%%%%%%%%%----------------
            if t < tDarcyOn

                % Zero Darcy velocities
                u(:) = 0;
                w(:) = 0;

                dt_dff = (h/2)^2 / max(kT(:)./rho(:)./cP(:) + eps);                                              
                dt     = CFL * dt_dff;               % time step [s]

                % Temperature update (no advection)
                dTdt = diffusion(T,kT,h,ix3,iz3,qbot)./(rho.*cP) ...
                     + Qr./(rho.*cP);
        
                T = T + dTdt * dt;
                
                % Enforce surface / air BC
                T(air) = Ttop;
                %T(1,:) = Ttop;


            %**********
            % After x years turn Darcy on
            %************
            else
                % initiate small perturbation when Darcy starts
                if t >= tDarcyOn && t < tDarcyOn + dt
                    T(~air) = T(~air) + 0.1*randn(nnz(~air),1);
                end

                % update pseudo-transient iterative step
                dtau = (h/2)^2 ./ max(KD(:) + eps);

                dTdto = dTdt;
                itP     = 0;
                res_rms = 1;
                dp      = 0*p;

                while res_rms >= tolP
    
                    % Temperature dependant density
                    rhoW = rhoW0 .* (1 - alphaT.*(T - Tref));
                    Drho = rhoW - mean(rhoW,2);
    
                    % Get Darcy flux, residuals, close boundaries
                    [u,w,resP] = darcy_flux(p, Drho, g, KD, h, ix3, iz3);
                    w(2:end,:) = (~air).*w(2:end,:); 
    
                    % pressure update 
                    dp = - alpha * dtau * resP + beta*dp;
                    p = p + dp;
    
                    % Root mean square error of residuals
                    den = max(rms(p(:)), 1e-10);
                    res_rms = rms(dtau*resP(:))/ den;
    
                    % print model progress
                    if ~mod(itP,nop)
                        fprintf('\n\n')
                        fprintf('*************************************************************\n');
                        fprintf('*****  GEOTHERMAL MODEL | %s  **************\n',datetime('now'));
                        fprintf('*************************************************************\n');
                        % Print time step information
                        fprintf(1,'*****  step %d;  dt = %4.4e; dtau = %4.4e;  time = %4.4e [yr]\n\n',k,dt/yr,dtau/yr,t/yr)
                        fprintf('\n   run ID: %s \n\n',runID);
                        fprintf(1,'  ---  iter %d;  res norm = %4.4e \n',itP,res_rms);
                        fprintf("T min/max: %g %g\n", min(T(~air),[],'all'), max(T(~air),[],'all'));
                       % fprintf("P min/max: %g %g\n", min(p(~air),[],'all'), max(p(~air),[],'all'));
                        q_est = -kT(end,:) .* (T(end,:) - T(end-1,:)) / h;
                        disp([min(q_est), mean(q_est), max(q_est)])
                    end

                    if itP >= 1e6
                        break
                    end
    
                    itP = itP + 1;
            end

                
            % update Temperature
            dTdt = diffusion(T,kT,h,ix3,iz3,qbot)./(rho.*cP)  ...
                 + advection(T,u,w,h,ix5,iz5,ADVN) ...
                 + Qr./(rho.*cP);
            
            % set time step size (Courant–Friedrichs–Lewy condition)
            dt_adv = (h/2)   / (max(abs(u(:))) + max(abs(w(:))) + eps); 
            dt_dff = (h/2)^2 / max(kT(:)./rho(:)./cP(:) + eps);                                              
            dt     = CFL * min(dt_adv,dt_dff);               % time step [s]
            
            resT = T - To - (dTdto+dTdt)/2 * dt; 
            
            % fix air layer temperature 
            resT(air) = 0;
            T(air) = Ttop;
            %T(1,:)   = Ttop;

            T = T - resT/3;
                end 

            % plot model progress
            if ~mod(k,nop)
                figure(3); clf
                Tempfig(xc,zc,T,t/yr); fig3=figure(3);
                figure(4); clf
                isotherm2D(xc,zc,T,t/yr); fig4=figure(4);
                figure(5); clf
                plotdrill(xc,zc,T,air,x_dh,zdrill,Tdrill,t,yr); fig5=figure(5);
                figure(6); clf
                isothermplot(xc,zc,T,air,x_pds,t,yr); fig6=figure(6);
                pause(0.1);
            end

            % set file name to runID_figtype_(frame#).mat and place in ../out/runID/
            % save figure handle 'fig1' every '4*nop' time steps
            if ~mod(k,20*nop)
             fig_3      = ['../out/',runID,'/',runID,'_','Tempfig',     '_',num2str(floor(k/(20*nop)))];
             fig_4   = ['../out/',runID,'/',runID,'_','isotherm2D',  '_',num2str(floor(k/(20*nop)))];
             fig_5    = ['../out/',runID,'/',runID,'_','plotdrill',   '_',num2str(floor(k/(20*nop)))];
             fig_6 = ['../out/',runID,'/',runID,'_','isothermplot','_',num2str(floor(k/(20*nop)))];
             print(fig3,fig_3,'-dpng','-r300','-image');
             print(fig4,fig_4,'-dpng','-r300','-image');
             print(fig5,fig_5,'-dpng','-r300','-image');
             print(fig6,fig_6,'-dpng','-r300','-image');
            end

    end

    % increment time
    t = t + dt;
end


switch MODE
    case 'VERIFY'
    % *****  calculate and display numerical error norm
    Ta  = analytical(T0,dT,sgm0,k0,u0,w0,Xc,Zc,D,W,t);
    Err = rms(T - Ta, 'all') / rms(Ta, 'all');
    
    disp(' ');
    disp(['Advection scheme: ',ADVN]);
    disp(['Time integration scheme: ',TINT]);
    disp(['Numerical error norm = ',num2str(Err)]);
    disp(' ');
end                         



%**************************************************************************
%*****  Utility Functions  ************************************************





%**************************************************************************
%*****  Function to calculate diffusion rate

function dfdt = diffusion(f,k,h,ix,iz,qbot)

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

% impose heat flux at boundary
qz(end,:) = -qbot;

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

% output variable
% dfdt: advection rate of scalar field


%****** horizontal advection (x direction) ******

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
q_ip_pos = u_pos(:,2:end).*f_ip_pos;          % flux on right face i+1/2 
q_im_pos = u_pos(:,1:end-1).*f_im_pos;        % flux on left face  i-1/2

% negative velocity
q_ip_neg = u_neg(:,2:end).*f_ip_neg;          % flux on right face i+1/2 
q_im_neg = u_neg(:,1:end-1).*f_im_neg;        % flux on left face  i-1/2

% advection flux balance for rate of change
div_qx_pos = (q_ip_pos - q_im_pos)/h;  % positive velocity                  
div_qx_neg = (q_ip_neg - q_im_neg)/h;  % negative velocity

div_qx     = div_qx_pos + div_qx_neg;     % combined flux x direction


%****** vertical advection (z direction) ********

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
dfdt  = - div_q;                        % advection rate
end


% Function for Darcy flux**********************************
function [u,w,res] = darcy_flux(p, Drho, g, KD, h, ix, iz)
 
% calculate Darcy flux coefficient at cell faces using harmonic averaging
KD1 = KD(iz(1:end-1),:);
KD2 = KD(iz(2:end),:);
KDz = 2*KD1.*KD2 ./ (KD1 + KD2 + eps);

KD1 = KD(:,ix(1:end-1));
KD2 = KD(:,ix(2:end));
KDx = 2*KD1.*KD2 ./ (KD1 + KD2 + eps);

% calculate density coefficient at cell faces
Drhoz = (Drho(iz(1:end-1),:)+Drho(iz(2:end),:))/2;

% calculate diffusive flux of scalar field f
w = - KDz .* diff(p(iz,:),1,1)/h + KDz .* Drhoz*g;
u = - KDx .* diff(p(:,ix),1,2)/h;

% Closed flow boundaries
w(1,:) = 0; w(end,:) = 0;
u(:,1) = 0; u(:,end) = 0;

% calculate Darcy flux balance for rate of change
res = diff(w,1,1)/h + diff(u,1,2)/h;
end



%****function to calculate the analytical temperature at time t 
function Ta = analytical(f0,df,sgm0,k0,u0,w0,Xc,Zc,D,W,t)

sgmt = sqrt(sgm0^2 + 2*k0*t);
% sum each of the 9 combinations
Ta   = f0 + df*(sgm0.^2/sgmt.^2)*(exp(-((Xc-(W/2    )- u0*t).^2 + (Zc-(D/2    ) - w0*t).^2)./ (4*sgmt^2)) ...
                                + exp(-((Xc-(W/2 + W)- u0*t).^2 + (Zc-(D/2    ) - w0*t).^2)./ (4*sgmt^2)) ...
                                + exp(-((Xc-(W/2 - W)- u0*t).^2 + (Zc-(D/2    ) - w0*t).^2)./ (4*sgmt^2)) ...
                                + exp(-((Xc-(W/2    )- u0*t).^2 + (Zc-(D/2 + D) - w0*t).^2)./ (4*sgmt^2)) ...
                                + exp(-((Xc-(W/2    )- u0*t).^2 + (Zc-(D/2 - D) - w0*t).^2)./ (4*sgmt^2)) ...
                                + exp(-((Xc-(W/2 + W)- u0*t).^2 + (Zc-(D/2 + D) - w0*t).^2)./ (4*sgmt^2)) ...
                                + exp(-((Xc-(W/2 + W)- u0*t).^2 + (Zc-(D/2 - D) - w0*t).^2)./ (4*sgmt^2)) ...
                                + exp(-((Xc-(W/2 - W)- u0*t).^2 + (Zc-(D/2 + D) - w0*t).^2)./ (4*sgmt^2)) ...
                                + exp(-((Xc-(W/2 - W)- u0*t).^2 + (Zc-(D/2 - D) - w0*t).^2)./ (4*sgmt^2)));
end







%*****  Function to make Temperature output figure 
function Tempfig(x,z,T,t)

imagesc(x,z,T); axis equal tight; colorbar
ylabel('z [m]','FontSize',15)
title(['Temp [C]; time = ',num2str(t),' [yr]'],'FontSize',17)
end


%*****  Function to make verification figure (verify against analytical solution)
function makefigverify(x,z,T,Ta,t) 
subplot(2,2,1);
imagesc(x,z,T); axis equal tight; colorbar
ylabel('z [m]','FontSize',15)
title(['Temperature [C]; time = ',num2str(t),' [yr]'],'FontSize',17)

subplot(2,2,2)
imagesc(x,z,T-Ta); axis equal tight; colorbar
xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('Error: T - Ta','FontSize',17)
drawnow;
end

%*****  Function to make output figure of initial values
function initialfig(x,z,T,KD,Qr,kT,t)
subplot(2,2,1);
imagesc(x,z,T); axis equal tight; colorbar
ylabel('z [m]','FontSize',15)
title(['Temperature [C]; time = ',num2str(t),' [yr]'],'FontSize',17)
subplot(2,2,2)
imagesc(x,z,KD); axis equal tight; colorbar
xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('Darcy mobility [m2/Pas]','FontSize',17)
subplot(2,2,3)
imagesc(x,z, Qr); axis equal tight; colorbar
xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('Radiogenic heating [W/m^3]','FontSize',17)
subplot(2,2,4)
imagesc(x,z,kT); axis equal tight; colorbar
xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('Thermal conductivity [W/m/K]','FontSize',17)
drawnow;
end

%***** Function to plot the model profile at a given x coordinate
% Granite Drillhole approx x = 5000 plus minus 100
% Proposed drillsite approx x = 11200

function plotdrill(xc,zc,T,air,x_dh,zdrill,Tdrill,t,yr)

% Find nearest model column to drill location
[~, ix] = min(abs(xc - x_dh));
    
% Extract full column
Tcol = T(:,ix);

% Find first non-air cell in this column
k0 = find(air(:,ix)==0, 1, 'first');   % first rock index

% Trim to rock-only and shift depth so surface = 0
Tmodel = Tcol(k0:end);
zrock  = zc(k0:end) - zc(k0);          % now starts at 0 at ground

% Create figure
plot(Tmodel, zrock, 'b-', 'LineWidth', 2); hold on
plot(Tdrill, zdrill, 'ro', 'MarkerFaceColor','r');

set(gca,'YDir','reverse')
xlabel('Temperature [°C]')
ylabel('Depth [m]')
title(sprintf('Drillhole Comparison at x ≈ %.0f m | t = %.1f yr', ...
              xc(ix), t/yr))
legend('Model (rock only)','Drillhole data')

% Interpolate model onto drill depths
Tmod_at_drill = interp1(zrock, Tmodel, zdrill, 'linear', 'extrap');

% Residuals
res  = Tmod_at_drill - Tdrill;

% Root mean squared error
RMSE = sqrt(mean(res.^2));

% Show RMSE on plot
txt = sprintf('RMSE = %.2f °C', RMSE);
text(0.04, 0.6, txt, 'Units','normalized', 'VerticalAlignment','bottom', ...
     'BackgroundColor','w', 'EdgeColor', 'k');
grid on
drawnow
end






% function to create isotherm depths at drillsite (with copilot assistance)
function isothermplot(xc,zc,T,air,x_pds,t,yr)

% Find nearest model column to proposed drill site
[~, ix] = min(abs(xc - x_pds));

Tcol = T(:,ix);
aircol = air(:,ix);

% surface = first rock cell (first non-air)
kSurf = find(~aircol, 1, 'first');

zRel = zc - zc(kSurf);          % depth below ground surface
Trock = Tcol(kSurf:end);
zRock = zRel(kSurf:end);

levels = [50 70 100 120];
zIso = nan(size(levels));

for k = 1:numel(levels)
    L = levels(k);
    j = find(Trock(1:end-1) < L & Trock(2:end) >= L, 1, 'first');
    if ~isempty(j)
        zIso(k) = interp1(Trock(j:j+1), zRock(j:j+1), L);
    elseif Trock(1) >= L
        zIso(k) = 0;  % already above at surface
    end
end

plot(Trock, zRock, 'b-', 'LineWidth', 2); hold on
set(gca,'YDir','reverse')
grid on

% only draw ylines if the isotherm exists
for k = 1:numel(levels)
    if ~isnan(zIso(k))
        yline(zIso(k), '--r', sprintf('%d°C @ %.0f m', levels(k), zIso(k)));
    end
end

xlabel('Temperature [°C]')
ylabel('Depth below ground [m]')
title(sprintf('Proposed drillsite at x ≈ %.0f m | t = %.1f yr', xc(ix), t/yr))
legend('Model','Location','southeast')
drawnow
end


% function to create 2D isotherm depths across domain 
function isotherm2D(x, z, T, t)
imagesc(x,z,T);
axis equal tight
set(gca,'YDir','reverse')
colorbar
hold on

% --- Add isotherm contour lines ---
levels = [50 70 100 120];   
[C,hc] = contour(x,z,T,levels,'k','LineWidth',1.5);

clabel(C,hc,'FontSize',10,'Color','k')

xlabel('x [m]')
ylabel('z [m]')
title(['Temp [C]; time = ',num2str(t),' [yr]'])

end