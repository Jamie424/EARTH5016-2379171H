%*****  1D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT  *******************


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
        ix3 = [      Nx, 1:Nx, 1   ];
        ix5 = [Nx-1, Nx, 1:Nx, 1, 2];
        iz3 = [      Nz, 1:Nz, 1   ];
        iz5 = [Nz-1, Nz, 1:Nz, 1, 2]; 
    case 'insulating'
        % example non-periodic indexing for N=4 
        % 3-point stencil            |-- i-1 --|-- i --|-- i+1 --|
        % 5-point stencil  |-- i-2 --|-- i-1 --|-- i --|-- i+1 --|-- i+2 --|
        ix3 = [     1, 1:Nx, Nx    ]; 
        ix5 = [1,   1, 1:Nx, Nx, Nx];  
        iz3 = [     1, 1:Nz, Nz    ]; 
        iz5 = [1,   1, 1:Nz, Nz, Nz];
end

% set initial coefficient fields
kT  = kT0  .* ones(Nz,Nx);         % matrix of initial thermal conductivity [W/m/K]
cP  = cP0  .* ones(Nz,Nx);         % heat capacity [J/kg/K]
rho = rho0 .* ones(Nz,Nx);         % density [kg/m3]
Qr  = Qr0  .* ones(Nz,Nx);         % heat productivity [W/m3]

% set initial velocity field
w = w0 .* ones(Nz+1,Nx);          % advection x-speed [m/s]
u = u0 .* ones(Nz,Nx+1);          % advection z-speed [m/s]

% set time step size (Courant–Friedrichs–Lewy condition)
dt_adv = (h/2)   / max(abs(u(:)) + abs(w(:)) + eps); 
dt_dff = (h/2)^2 / max(kT(:)./rho(:)./cP(:)+eps);
CFL = 0.5;                                               
dt     = CFL * min(dt_adv,dt_dff); % time step [s]

% set initial temperature field
T = analytical(T0,dT,sgm0,kT0./rho0./cP0,u0,w0,Xc,Zc,D,W,t);   % analytical solution in 2D
%T   = T0 + dT*exp(-(xc-W/2).^2./(2*sgm0^2));      
Tin = T;                                          % store initial condition for plotting
Ta  = T;                                          % initialise analytical solution

% Initialise time count variables
t = 0;  % initial time [s]
k = 0;  % initial time step count

% initialise output figure with initial condition
figure(1); clf
makefig(xc,T,Tin,Ta,0);


%********** Implicit scheme setup

% set temporal coefficient matrix
At = speye(N,N) * 1/dt; % N x N sparse diagonal matrix with 1/dt on diagonal

% use mapping array to add indices to index lists for spatial coefficients 
i = []; j = []; % initialise i,j as empty lists 
i = [i, ind3(2:end-1)]; j = [j, ind3(2:end-1)]; % centre stencil node i 
i = [i, ind3(2:end-1)]; j = [j, ind3(1:end-2)]; % left stencil node i-1 
i = [i, ind3(2:end-1)]; j = [j, ind3(3:end  )]; % right stencil node i+1 

% add coefficient values to value list 
a = []; % initialise a as empty list 
a = [a,  ( +2*k0/dx^2)        * ones(1,N)];       % centre stencil node i 
a = [a, -(u0/(2*dx)- k0/dx^2) * ones(1,N)];       % left stencil node i-1
a = [a,  (u0/(2*dx)- k0/dx^2) * ones(1,N)];       % right stencil node i+1

disp([numel(i), numel(j), numel(a)])  % length of i,j and a must be equal to pair the three together in Ax

Ax = sparse(i,j,a,N,N);   % place values a at positions (i,j) in NxN matrix
A_BE = At + Ax;           % coefficient matrix for BE1 scheme
A_CN = At + Ax/2;         % coefficient matrix for CN2 scheme

%*****  Solve Model Equations
while t <= tend

    % increment time and step count
    t = t+dt;
    k = k+1;
    
    
    switch SCHEME
        case 'explicit' 
            switch TINT % select time integration scheme
                case 'FE1'  % 1st-order Forward Euler time integration scheme
                    % get rate of change
                    dTdt = ((diffusion(T,k,h,ix,iz3) + Qr0)./(rho0*cp0))  ...
                           - advection(T,u0,w0,h,iz5,ADVN);
        
                case 'RK2'  % 2nd-order Runge-Kutta time integration scheme
                    dTdt_half = diffusion(T               ,k,h,ix,iz) + Qr0 ./(rho0*cp0) ...
                              + advection(T,u0,dx,ind5,ADVN);
                    dTdt      = diffusion(T+dTdt_half*dt/2,k0,dx,iz3) ...
                              + advection(T+dTdt_half*dt/2,u0,dx,iz5,ADVN);
            end
            % update temperature
            T = T + dTdt * dt;

        case 'implicit'
            switch TINT % select time integration scheme
                case 'BE1' % implicit: 1st-order Backward Euler (BE1) 
                    b = At*T.';            % prepare forcing vector 
                    T = (A_BE \ b).';      % solve linear system of equations

                case 'CN2' % Semi-implicit: 2nd-order Crank-Nicolson (CN2) 
                    b = At*T.'- Ax*T.'/2;   % get forcing vector
                    T = (A_CN \ b).';       % solve linear system of equations 
            end    
    end


    % get analytical solution at time t
    analytical(T0,df,sgm0,k0,u0,w0,Xc,Zc,D,W,t)

    % plot model progress
    if ~mod(k,nop)
        makefig(xc,T,Tin,Ta,t/yr);
        pause(0.1);
    end




%*****  calculate and display numerical error norm
Err = rms(T - Ta, 'all') / rms(Ta, 'all');

disp(' ');
disp(['Advection scheme: ',ADVN]);
disp(['Time integration scheme: ',TINT]);
disp(['Numerical error = ',num2str(Err)]);
disp(' ');

end

%*****  Utility Functions  ************************************************

%*****  Function to make output figure

function makefig(x,T,Tin,Ta,t)

subplot(2,1,1)
plot(x,Tin,'k:',x,T,'r-',x,Ta,'k--','LineWidth',1.5); axis tight; box on;

ylabel('T [C]','FontSize',15)
title(['Temperature at time = ',num2str(t,4),' yr'],'FontSize',18)

subplot(2,1,2)
plot(x,(T-Ta)./rms(Ta,'all'),'r-',x,0*Ta,'k-','LineWidth',1.5); axis tight; box on;

xlabel('x [m]','FontSize',15)
ylabel('E [1]','FontSize',15)
title(['Numerical Error at time = ',num2str(t,4),' yr'],'FontSize',18)

drawnow;

end


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


%*****  Function to calculate advection rate

function dfdt = advection(f,u,w,h,ix,iz,ADVN)

% input arguments
% f:    advected scalar field
% u:    advection velocity x direction 
% w:    advection velocity z direction
% h:    grid spacing
% ix:  ghosted index list x direction
% iz:  ghosted index list z direction
% ADVN: advection scheme

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
q_ip_pos = u_pos.*f_ip_pos;        % flux on right face i+1/2 
q_im_pos = u_pos.*f_im_pos;        % flux on left face  i-1/2

% negative velocity
q_ip_neg = u_neg.*f_ip_neg;        % flux on right face i+1/2 
q_im_neg = u_neg.*f_im_neg;        % flux on left face  i-1/2

% advection flux balance for rate of change
div_qx_pos = (q_ip_pos - q_im_pos)/h;  % positive velocity                  
div_qx_neg = (q_ip_neg - q_im_neg)/h;  % negative velocity

div_qx     = div_q_pos + div_q_neg;     % combined flux x direction

div_q = div_qx + div_qz     % xflux + z flux
dfdt  = - div_q;            % advection rate



%****** horizontal advection (x direction)


% split the velocities into positive and negative
w_pos = max(0,w);    % positive velocity (to the right)
w_neg = min(0,w);    % negative velocity (to the left)


f_jmm = f(iz(1:end-4),:);   % i - 2
f_jm  = f(iz(2:end-3),:);   % i - 1
f_jc  = f(iz(3:end-2),:);   % i
f_jp  = f(iz(4:end-1),:);   % i + 1
f_jpp = f(iz(5:end  ),:);   % i + 2




end

%****function to calculate the analytical temperature at time t
function Ta = analytical(f0,df,sgm0,k0,u0,w0,Xc,Zc,D,W,t)

sgmt = sqrt(sgm0^2 + 2*k0*t);
% sum each of the 9 combinations
Ta   = f0 + dT*(sgm0/sgmt)*(exp(-(Xc-(W/2    )- u0*t).^2 + (Zc-(D/2    ) - w0*t).^2)./ (2*sgmt^2)) ...
                          + exp(-(Xc-(W/2 + W)- u0*t).^2 + (Zc-(D/2    ) - w0*t).^2)./ (2*sgmt^2)) ...
                          + exp(-(Xc-(W/2 - W)- u0*t).^2 + (Zc-(D/2    ) - w0*t).^2)./ (2*sgmt^2)) ...
                          + exp(-(Xc-(W/2    )- w0*t).^2 + (Zc-(D/2 + D) - w0*t).^2)./ (2*sgmt^2)) ...
                          + exp(-(Xc-(W/2    )- u0*t).^2 + (Zc-(D/2 - D) - w0*t).^2)./ (2*sgmt^2)) ...
                          + exp(-(Xc-(W/2 + W)- u0*t).^2 + (Zc-(D/2 + D) - w0*t).^2)./ (2*sgmt^2)) ...
                          + exp(-(Xc-(W/2 + W)- u0*t).^2 + (Zc-(D/2 - D) - w0*t).^2)./ (2*sgmt^2)) ...
                          + exp(-(Xc-(W/2 - W)- u0*t).^2 + (Zc-(D/2 + D) - w0*t).^2)./ (2*sgmt^2))...
                          + exp(-(Xc-(W/2 - W)- u0*t).^2 + (Zc-(D/2 - D) - w0*t).^2)./ (2*sgmt^2)))

end
%test


% get analytical function to work

% 