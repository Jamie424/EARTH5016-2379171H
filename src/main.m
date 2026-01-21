%*****  1D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT  *******************


%*****  Initialise Model Setup

% create x-coordinate vectors
xc = dx/2:dx:W-dx/2;    % coordinate vector for cell centre positions [m]
xf = 0:dx:W;            % coordinate vectore for cell face positions [m]

% set time step size
dt_adv = (dx/2)   / u0;
dt_dff = (dx/2)^2 / k0;
dt     = CFL * min(dt_adv,dt_dff); % time step [s]

% set up ghosted index lists for boundary conditions
switch BC
    case 'periodic'
        ind3 = [     N, 1:N, 1   ];  % 3-point stencil            |-- i-1 --|-- i --|-- i+1 --|
        ind5 = [N-1, N, 1:N, 1 ,2];  % 5-point stencil  |-- i-2 --|-- i-1 --|-- i --|-- i+1 --|-- i+2 --|
    case 'insulating'
        % example non-periodic indexing for N=4 
        ind3 = [     1, 1:N, N   ];  % 3-point stencil            |-- i-1 --|-- i --|-- i+1 --|
        ind5 = [1,   1, 1:N, N, N];  % 5-point stencil  |-- i-2 --|-- i-1 --|-- i --|-- i+1 --|-- i+2 --|
end

% set initial condition for temperature at cell centres
T   = T0 + dT*exp(-(xc-W/2).^2./(2*sgm0^2));      % initialise T array at Tr
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
a = [a, -(u0/(2*dx)- k0/dx^2) * ones(1,N)];       % left stencil node i-1   Note: vector of 1s length N-2 as ind3(2:end-1) is length N-2
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
                    dTdt = diffusion(T,k0,dx,ind3) ...
                         + advection(T,u0,dx,ind5,ADVN);
        
                case 'RK2'  % 2nd-order Runge-Kutta time integration scheme
                    dTdt_half = diffusion(T               ,k0,dx,ind3) ...
                              + advection(T               ,u0,dx,ind5,ADVN);
                    dTdt      = diffusion(T+dTdt_half*dt/2,k0,dx,ind3) ...
                              + advection(T+dTdt_half*dt/2,u0,dx,ind5,ADVN);
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
    sgmt = sqrt(sgm0^2 + 2*k0*t);
    Ta   = T0 + dT*(sgm0/sgmt)*exp(-(xc -(W/2)    - u0*t).^2 ./ (2*sgmt^2)) ...
              + dT*(sgm0/sgmt)*exp(-(xc -(W/2 + W)- u0*t).^2 ./ (2*sgmt^2)) ...
              + dT*(sgm0/sgmt)*exp(-(xc -(W/2 - W)- u0*t).^2 ./ (2*sgmt^2));

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

function dfdt = diffusion(f,k,dx,ind)

% input arguments
% f:    diffusing scalar field
% k:    diffusion coefficient
% dx:   grid spacing
% ind:  ghosted index list

% output variables
% dfdt: diffusion rate of scalar field f

% calculate diffusive flux of scalar field f
q = - k * diff(f(ind)) / dx;

% calculate diffusion flux balance for rate of change
dfdt = - diff(q) / dx;

% Note: start of function defines the output dfdt
end


%*****  Function to calculate advection rate

function dfdt = advection(f,u,dx,ind,ADVN)

% input arguments
% f:    advected scalar field
% u:    advection velocity
% dx:   grid spacing
% ind:  ghosted index list
% ADVN: advection scheme

% output variables
% dfdt: advection rate of scalar field

% split the velocities into positive and negative
u_pos = max(0,u);    % positive velocity (to the right)
u_neg = min(0,u);    % negative velocity (to the left)

% get values on stencil nodes
f_imm  = f(ind(1:end-4));  % i-2
f_im   = f(ind(2:end-3));  % i-1   
f_ic   = f(ind(3:end-2));  % i
f_ip   = f(ind(4:end-1));  % i+1
f_ipp  = f(ind(5:end  ));  % i+2

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
div_q_pos = (q_ip_pos - q_im_pos)/dx;  % positive velocity                  
div_q_neg = (q_ip_neg - q_im_neg)/dx;  % negative velocity

div_q     = div_q_pos + div_q_neg;     % combined
dfdt      = - div_q;                   % advection rate

end



%test