function F_int = compute_fluxes_HLL(U)

% This function computes all interface numerical fluxes, using 
% the HLL scheme.

phys_const; % Load physical constants

N_int = size(U,2)-1; % aka N_int = N_cells - 1

% Loop on interfaces
for ii = 1:N_int

  % ID CELLS
  ID_C_L = ii;
  ID_C_R = ii+1;
  
  % Solution at the Left and Right of the interface
  UL = U(:, ID_C_L);
  UR = U(:, ID_C_R);
  
  % Reconstruct solution (aka primitive variables)
  rhoL = UL(1);
  uL   = UL(2)/rhoL;
  PL   = (gam-1)*(UL(3) - rhoL*uL*uL/2);
  
  rhoR = UR(1);
  uR   = UR(2)/rhoR;
  PR   = (gam-1)*(UR(3) - rhoR*uR*uR/2);
  
  % Build flux function (from the Euler equations)
  FL = [0;0;0];
  FL(1) = rhoL*uL;
  FL(2) = rhoL*uL*uL + PL;
  FL(3) = rhoL*uL^3/2 + PL/(gam-1)*uL + PL*uL;
  
  FR = [0;0;0];
  FR(1) = rhoR*uR;
  FR(2) = rhoR*uR*uR + PR;
  FR(3) = rhoR*uR^3/2 + PR/(gam-1)*uR + PR*uR;
  
  % Compute wave speeds
  aL = sqrt(gam*PL/rhoL); % Speed of sound in the left cell
  aR = sqrt(gam*PR/rhoR); % Speed of sound in the right cell
  
  sMinL = uL - aL; % Minimum wave speed (Euler's eq eigenvalues)
  sMaxR = uR + aR; % Maximum wave speed (Euler's eq eigenvalues)
 
  % HLL fluxes 
  if (sMinL >= 0)
    F_int(:,ii) = FL;
  elseif (sMaxR <= 0)
    F_int(:,ii) = FR;
  else
    F_int(:,ii) = (sMinL*sMaxR*(UR - UL) + sMaxR*FL - sMinL*FR)/(sMaxR - sMinL);
  end

end
