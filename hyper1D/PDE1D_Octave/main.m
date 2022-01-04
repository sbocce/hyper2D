close all
clear
clc

% For Octave
try
  page_screen_output(0);
end

% Load phys stuff
phys_const;

% =============== Time discretization ====================
dt     = 2e-6;     % [s] time step 
N_time = 200;      % Number of time steps
t_vect = [0:dt:dt*N_time]; % Array of times

% =============== Create mesh ====================
x_min   = 0.0;         % [m] Domain size
x_max   = 1.0;         % [m] Domain size
N_cells = 500;         % Number of cells, including ghost cells
N_int   = N_cells - 1; % Number of interfaces

Lcell  = (x_max - x_min)/(N_cells-2); % [m] size of a single cell
x_cc   = linspace(x_min - Lcell/2, x_max + Lcell/2, N_cells);
x_int  = (x_cc(2:end) + x_cc(1:end-1))/2;

% ============== Solution vector ==================

% The solution vector U contains 3 rows (one for the density, one for the momentum
% and one for the energy) and one column for each grid cell 

U = zeros(3, N_cells); 

% +++++ Initial condition - initialize as the Sod shocktube problem

% Left state
rhoL = 4.0; % [kg/m3] density
uL   = 0.0; % [m/s] velocity
PL   = 400000; % [Pa] pressure

U(1,1:floor(N_cells/2)) = rhoL;
U(2,1:floor(N_cells/2)) = rhoL*uL;
U(3,1:floor(N_cells/2)) = rhoL*uL^2/2 + PL/(gam-1);

% Right state
rhoR = 1.0; % [kg/m3] density 
uR   = 0.0; % [m/s] velocity
PR   = 100000; % [Pa] pressure

U(1,floor(N_cells/2):end) = rhoR;
U(2,floor(N_cells/2):end) = rhoR*uR;
U(3,floor(N_cells/2):end) = rhoR*uR^2/2 + PR/(gam-1);

% =============== SOLVE IN TIME ===============
SolFig = figure();

for tID = 1:N_time
  
  fprintf('Solving timestep %d of %d\n', tID, N_time);
    
  % Compute interface numerical fluxes from solution  
  F_int = compute_fluxes_HLL(U); % This is a matrix, with N_int columns

  % Update solution (Forward Euler)
  % The solution in a cell is the solution itself PLUS the flux from the left interface
  % MINUS the flux from the right interface
  U(:,2:end-1) = U(:, 2:end-1) + dt/Lcell*(F_int(:, 1:end-1) - F_int(:, 2:end));
  
  % Plot the solution every Nplot steps
  Nplot = 10;
  if (mod(tID,Nplot)==0) 
    subplot(3,1,1)
    plot(x_cc, U(1,:),'b','linewidth',2);
    xlabel('Position [m]')
    ylabel('Density [kg/m3]')
    title(['Time: ', num2str(t_vect(tID)), ' s']);
  
    subplot(3,1,2) 
    plot(x_cc, U(2,:),'b','linewidth',2);
    xlabel('Position [m]')
    ylabel('Momentum [kg/m3*m/s]')
  
    subplot(3,1,3)
    plot(x_cc, U(3,:),'b','linewidth',2);
    xlabel('Position [m]')
    ylabel('Energy [J/m3]')
    pause(0.001)
  end
  
end
