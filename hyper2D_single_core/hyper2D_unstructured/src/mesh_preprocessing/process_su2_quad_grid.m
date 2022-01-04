%%% 
% It recognizes "inlet" "outlet" "wall" and "sym" (symmetry) BCs. Use EXACLTY these markers in gmsh.
% 

close all
clear
clc

try
  page_screen_output(0);
end

% Map inlet names to IDs
BC_inlet_ID  = -1;
BC_outlet_ID = -2;
BC_wall_ID   = -3;
BC_sym_ID    = -4; 

% ------
% Init some stuff (working variables)
inlet_lines  = [-1,-1];
outlet_lines = [-1,-1];
wall_lines   = [-1,-1];
sym_lines    = [-1,-1];

% ########### Load mesh
fid = fopen('cyl_mesh_finer.su2','r');
% fid = fopen('cyl_finefine.su2','r');
% fid = fopen('mesh_deform.su2','r');
% fid = fopen('mesh_ref.su2','r');
% fid = fopen('mesh_quad_test.su2','r');

% Loop on file lines 
while( (str = fgetl(fid)) ~= -1) % Read until the end

  fprintf('%s\n', str)

  % +++++++ Read elements ++++++++++
   
  if strcmp(str(1:6), 'NELEM=')

    [dummy, Nelem] = sscanf(str, '%s %d', 'C')

    % Create array 
    ele_list_quad = zeros(Nelem,4); 

    for ele_ID = 1:Nelem

      str = fgetl(fid); % Read new line

      % syntax: 
      % eleType | point1 | point2 | point3 | point4 | eleID
      all_int = sscanf(str, '%d %d %d %d %d %d');

      if (all_int(1) == 9) % QUAD

        % SU2 counts starting from 0. Add 1
        p1 = all_int(2) + 1;
        p2 = all_int(3) + 1;
        p3 = all_int(4) + 1;
        p4 = all_int(5) + 1;

        ele_list_quad(ele_ID,:) = [p1, p2, p3, p4];

      else

        error('Only 2D quads supported!\n')

      end 

    end

  end % end read elements

  % ++++++++++ Read nodes ++++++++++

  if strcmp(str(1:6), 'NPOIN=')

    [dummy, Nnodes] = sscanf(str, '%s %d', 'C')

    % Create array of nodes
    nodes_list = zeros(Nnodes,2); % (x,y) of each node

    for node_ID = 1:Nnodes

      str = fgetl(fid); % Read new line

      all_float = sscanf(str, '%f %f %d');

      x_now = all_float(1);
      y_now = all_float(2);

      nodes_list(node_ID,:) = [x_now, y_now];

    end

  end % end read nodes

  % ++++++++ Read BCs +++++++++

  % ----------------- INLET BCs

  if strcmp(str, 'MARKER_TAG= inlet') 

      str = fgetl(fid); % Read new line

      [dummy, Ninlet] = sscanf(str, '%s %d', 'C');

      inlet_lines = zeros(Ninlet, 2);

      for bound_ID = 1:Ninlet

        str = fgetl(fid); % Read new line

        int_read = sscanf(str, '%d %d');

        inlet_lines(bound_ID,:) = [int_read(2), int_read(3)] + 1; % Remember: points are indexed starting from 0 in su2 mesh
      end 
  end % end read BCs

  % ----------------- OUTLET BCs
  if strcmp(str, 'MARKER_TAG= outlet') 

      str = fgetl(fid); % Read new line

      [dummy, Noutlet] = sscanf(str, '%s %d', 'C');

      outlet_lines = zeros(Noutlet, 2);

      for bound_ID = 1:Noutlet

        str = fgetl(fid); % Read new line

        int_read = sscanf(str, '%d %d');

        outlet_lines(bound_ID,:) = [int_read(2), int_read(3)] + 1; % Remember: points are indexed starting from 0 in su2 mesh
      end 
  end % end read BCs

  % ----------------- WALL BCs
  if strcmp(str, 'MARKER_TAG= wall') 

      str = fgetl(fid); % Read new line

      [dummy, Nwall] = sscanf(str, '%s %d', 'C');

      wall_lines = zeros(Nwall, 2); 

      for bound_ID = 1:Nwall

        str = fgetl(fid); % Read new line

        int_read = sscanf(str, '%d %d');

        wall_lines(bound_ID,:) = [int_read(2), int_read(3)] + 1; % Remember: points are indexed starting from 0 in su2 mesh
      end 
  end % end read BCs

  % ----------------- SYMMETRY BCs
  if strcmp(str, 'MARKER_TAG= sym') 

      str = fgetl(fid); % Read new line

      [dummy, Nsym] = sscanf(str, '%s %d', 'C');

      sym_lines = zeros(Nsym, 2); 

      for bound_ID = 1:Nsym

        str = fgetl(fid); % Read new line

        int_read = sscanf(str, '%d %d');

        sym_lines(bound_ID,:) = [int_read(2), int_read(3)] + 1; % Remember: points are indexed starting from 0 in su2 mesh
      end 
  end % end read BCs


end

xP = nodes_list(:,1);
yP = nodes_list(:,2);


% % PLOT ELEMENTS AND POINTS 
% figure
% plot(nodes_list(:,1), nodes_list(:,2), 'or', 'linewidth',2)
% axis equal
% for ii = 1:size(ele_list_quad)
% 
%   p1 = ele_list_quad(ii,1);
%   p2 = ele_list_quad(ii,2);
%   p3 = ele_list_quad(ii,3);
%   p4 = ele_list_quad(ii,4);
% 
% 
%   hold on
%   plot([nodes_list(p1,1), nodes_list(p2,1), nodes_list(p3,1), nodes_list(p4,1), nodes_list(p1,1)], [nodes_list(p1,2),nodes_list(p2,2), nodes_list(p3,2), nodes_list(p4,2), nodes_list(p1,2)], 'k', 'linewidth', 2)
% 
% %  pause(0.0001)
% end

% ++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++ WRITE NODES ++++++++++++++++
fid_nodes = fopen('nodes.hyp', 'w');

fprintf(fid_nodes, '%d\n', size(nodes_list,1));

for node_ID = 1:Nnodes
  fprintf(fid_nodes, '%e %e \n',nodes_list(node_ID,1), nodes_list(node_ID,2));
end

fclose(fid_nodes)

% ++++++++++++++++ Elaborate connectivity, compute element properties etc etc ++++++++++++++++++++

tic

fid_write = fopen('mesh.hyp', 'w');
fprintf(fid_write, '# ele_ID, A, xC, yC, L12, n12(1), n12(2), neigh_12, L23, n23(1), n23(2), neigh_23, L34, n34(1), n34(2), neigh_34, L41, n41(1), n41(2), neigh_41\n')
fprintf(fid_write, '%d\n', Nelem) % Write number of elements

ele_all = zeros(Nelem, 19);
for ele_ID = 1:Nelem

  fprintf('Processing element %d of %d\n', ele_ID, Nelem)

  p1 = ele_list_quad(ele_ID,1);
  p2 = ele_list_quad(ele_ID,2);
  p3 = ele_list_quad(ele_ID,3);
  p4 = ele_list_quad(ele_ID,4);

  % ===== FIND INTERFACES LENGTH =====

  L12 = sqrt( (xP(p2)-xP(p1)).^2 + (yP(p2)-yP(p1)).^2 );
  L23 = sqrt( (xP(p2)-xP(p3)).^2 + (yP(p2)-yP(p3)).^2 );
  L34 = sqrt( (xP(p3)-xP(p4)).^2 + (yP(p3)-yP(p4)).^2 );
  L41 = sqrt( (xP(p4)-xP(p1)).^2 + (yP(p4)-yP(p1)).^2 );

  % ===== COMPUTE AREA (use Heron's formula for the triangles, exploiting the interfaces length) ======
  % Assuming a convex quadrilateral, divide it in two triangles
  % 
  %       (4) O-------O (3)
  %           |\ ii  /
  %           | \   / 
  %           |i \ /
  %       (1) O---O (2)
  % 
  L24 = sqrt((xP(p4)-xP(p2)).^2 + (yP(p4)-yP(p2)).^2);
  s_i = (L12 + L24 + L41)/2;
  A_i = sqrt(s_i*(s_i-L12)*(s_i-L24)*(s_i-L41));

  s_ii = (L23 + L34 + L24)/2;
  A_ii = sqrt(s_ii*(s_ii-L23)*(s_ii-L34)*(s_ii-L24));

  A = A_i + A_ii;

  % ===== COMPUTE NORMALS =====
  v12 = [xP(p2) - xP(p1); yP(p2) - yP(p1); 0];
  v23 = [xP(p3) - xP(p2); yP(p3) - yP(p2); 0];
  v34 = [xP(p4) - xP(p3); yP(p4) - yP(p3); 0];
  v41 = [xP(p1) - xP(p4); yP(p1) - yP(p4); 0];

  n12 = cross(v12, [0;0;1]);
  n12 = n12/norm(n12);

  n23 = cross(v23, [0;0;1]);
  n23 = n23/norm(n23);

  n34 = cross(v34, [0;0;1]);
  n34 = n34/norm(n34);

  n41 = cross(v41, [0;0;1]);
  n41 = n41/norm(n41);

  % ===== FIND CENTROID =====

  xC = (xP(p1) + xP(p2) + xP(p3) + xP(p4))/4;
  yC = (yP(p1) + yP(p2) + yP(p3) + yP(p4))/4;

  % ===== FIND NEIGHBORS =====

  % the following arrays "side_ij" are made by at most two elements: 
  % 1) the cell itself
  % 2) if not a boundary element, a neighboring element
  % Note that the order may vary

  side_12 = find( ( ele_list_quad(:,1)==p1 | ele_list_quad(:,2)==p1 | ele_list_quad(:,3)==p1 | ele_list_quad(:,4)==p1 ) & ...
                  ( ele_list_quad(:,1)==p2 | ele_list_quad(:,2)==p2 | ele_list_quad(:,3)==p2 | ele_list_quad(:,4)==p2 ) );

  side_23 = find( ( ele_list_quad(:,1)==p2 | ele_list_quad(:,2)==p2 | ele_list_quad(:,3)==p2 | ele_list_quad(:,4)==p2 ) & ...
                  ( ele_list_quad(:,1)==p3 | ele_list_quad(:,2)==p3 | ele_list_quad(:,3)==p3 | ele_list_quad(:,4)==p3 ) );

  side_34 = find( ( ele_list_quad(:,1)==p3 | ele_list_quad(:,2)==p3 | ele_list_quad(:,3)==p3 | ele_list_quad(:,4)==p3 ) & ...
                  ( ele_list_quad(:,1)==p4 | ele_list_quad(:,2)==p4 | ele_list_quad(:,3)==p4 | ele_list_quad(:,4)==p4 ) );

  side_41 = find( ( ele_list_quad(:,1)==p4 | ele_list_quad(:,2)==p4 | ele_list_quad(:,3)==p4 | ele_list_quad(:,4)==p4 ) & ...
                  ( ele_list_quad(:,1)==p1 | ele_list_quad(:,2)==p1 | ele_list_quad(:,3)==p1 | ele_list_quad(:,4)==p1 ) );


  neigh_12 = sum(side_12) - ele_ID;  % By summing and removing ele_ID, we get the right neighbor ID
  neigh_23 = sum(side_23) - ele_ID;
  neigh_34 = sum(side_34) - ele_ID;
  neigh_41 = sum(side_41) - ele_ID;

  if neigh_12 == 0
    % TODO ASSIGN BC TYPE
    
    % ++++ Are the two points inlet boundaries?
    side_12 = find( ( inlet_lines(:,1) == p1 | inlet_lines(:,2) == p1 ) & ...
                    ( inlet_lines(:,1) == p2 | inlet_lines(:,2) == p2 ) );
    if (numel(side_12) == 1)  % Yes, found it!
      neigh_12 = BC_inlet_ID;
    end 

    % ++++ Are the two points outlet boundaries?
    side_12 = find( ( outlet_lines(:,1) == p1 | outlet_lines(:,2) == p1 ) & ...
                    ( outlet_lines(:,1) == p2 | outlet_lines(:,2) == p2 ) );
    if (numel(side_12) == 1)  % Yes, found it!
      neigh_12 = BC_outlet_ID;
    end 

    % ++++ Are the two points wall boundaries?
    side_12 = find( ( wall_lines(:,1) == p1 | wall_lines(:,2) == p1 ) & ...
                    ( wall_lines(:,1) == p2 | wall_lines(:,2) == p2 ) );
    if (numel(side_12) == 1)  % Yes, found it!
      neigh_12 = BC_wall_ID;
    end 

    % ++++ Are the two points symmetry boundaries?
    side_12 = find( ( sym_lines(:,1) == p1 | sym_lines(:,2) == p1 ) & ...
                    ( sym_lines(:,1) == p2 | sym_lines(:,2) == p2 ) );
    if (numel(side_12) == 1)  % Yes, found it!
      neigh_12 = BC_sym_ID;
    end 

  end
  
  % -------------------------
  if neigh_23 == 0 

    % ++++ Are the two points inlet boundaries?
    side_23 = find( ( inlet_lines(:,1) == p2 | inlet_lines(:,2) == p2 ) & ...
                    ( inlet_lines(:,1) == p3 | inlet_lines(:,2) == p3 ) );
    if (numel(side_23) == 1)  % Yes, found it!
      neigh_23 = BC_inlet_ID;
    end 

    % ++++ Are the two points outlet boundaries?
    side_23 = find( ( outlet_lines(:,1) == p2 | outlet_lines(:,2) == p2 ) & ...
                    ( outlet_lines(:,1) == p3 | outlet_lines(:,2) == p3 ) );
    if (numel(side_23) == 1)  % Yes, found it!
      neigh_23 = BC_outlet_ID;
    end 

    % ++++ Are the two points wall boundaries?
    side_23 = find( ( wall_lines(:,1) == p2 | wall_lines(:,2) == p2 ) & ...
                    ( wall_lines(:,1) == p3 | wall_lines(:,2) == p3 ) );
    if (numel(side_23) == 1)  % Yes, found it!
      neigh_23 = BC_wall_ID;
    end 

    % ++++ Are the two points symmetry boundaries?
    side_23 = find( ( sym_lines(:,1) == p2 | sym_lines(:,2) == p2 ) & ...
                    ( sym_lines(:,1) == p3 | sym_lines(:,2) == p3 ) );
    if (numel(side_23) == 1)  % Yes, found it!
      neigh_23 = BC_sym_ID;
    end 

  end

  % -------------------------
  if neigh_34 == 0 

    % ++++ Are the two points inlet boundaries?
    side_34 = find( ( inlet_lines(:,1) == p3 | inlet_lines(:,2) == p3 ) & ...
                    ( inlet_lines(:,1) == p4 | inlet_lines(:,2) == p4 ) );
    if (numel(side_34) == 1)  % Yes, found it!
      neigh_34 = BC_inlet_ID;
    end 

    % ++++ Are the two points outlet boundaries?
    side_34 = find( ( outlet_lines(:,1) == p3 | outlet_lines(:,2) == p3 ) & ...
                    ( outlet_lines(:,1) == p4 | outlet_lines(:,2) == p4 ) );
    if (numel(side_34) == 1)  % Yes, found it!
      neigh_34 = BC_outlet_ID;
    end 

    % ++++ Are the two points inlet boundaries?
    side_34 = find( ( wall_lines(:,1) == p3 | wall_lines(:,2) == p3 ) & ...
                    ( wall_lines(:,1) == p4 | wall_lines(:,2) == p4 ) );
    if (numel(side_34) == 1)  % Yes, found it!
      neigh_34 = BC_wall_ID;
    end 

    % ++++ Are the two points symmetry boundaries?
    side_34 = find( ( sym_lines(:,1) == p3 | sym_lines(:,2) == p3 ) & ...
                    ( sym_lines(:,1) == p4 | sym_lines(:,2) == p4 ) );
    if (numel(side_34) == 1)  % Yes, found it!
      neigh_34 = BC_sym_ID;
    end 

  end


  % -------------------------
  if neigh_41 == 0 

    % ++++ Are the two points inlet boundaries?
    side_41 = find( ( inlet_lines(:,1) == p4 | inlet_lines(:,2) == p4 ) & ...
                    ( inlet_lines(:,1) == p1 | inlet_lines(:,2) == p1 ) );
    if (numel(side_41) == 1)  % Yes, found it!
      neigh_41 = BC_inlet_ID;
    end 

    % ++++ Are the two points outlet boundaries?
    side_41 = find( ( outlet_lines(:,1) == p4 | outlet_lines(:,2) == p4 ) & ...
                    ( outlet_lines(:,1) == p1 | outlet_lines(:,2) == p1 ) );
    if (numel(side_41) == 1)  % Yes, found it!
      neigh_41 = BC_outlet_ID;
    end 

    % ++++ Are the two points wall boundaries?
    side_41 = find( ( wall_lines(:,1) == p4 | wall_lines(:,2) == p4 ) & ...
                    ( wall_lines(:,1) == p1 | wall_lines(:,2) == p1 ) );
    if (numel(side_41) == 1)  % Yes, found it!
      neigh_41 = BC_wall_ID;
    end 

    % ++++ Are the two points symmetry  boundaries?
    side_41 = find( ( sym_lines(:,1) == p4 | sym_lines(:,2) == p4 ) & ...
                    ( sym_lines(:,1) == p1 | sym_lines(:,2) == p1 ) );
    if (numel(side_41) == 1)  % Yes, found it!
      neigh_41 = BC_sym_ID;
    end 

  end

  % ======== ASSEMBLE ALL DATA 
  ele_all(ele_ID, :) = [A, xC, yC, L12, n12(1), n12(2), neigh_12, L23, n23(1), n23(2), neigh_23, L34, n34(1), n34(2), neigh_34, L41, n41(1), n41(2), neigh_41];

  fprintf(fid_write, '%d %e %e %e %e %e %e %d %e %e %e %d %e %e %e %d %e %e %e %d %d %d %d %d\n', ele_ID, A, xC, yC, L12, n12(1), n12(2), neigh_12, L23, n23(1), n23(2), neigh_23, L34, n34(1), n34(2), neigh_34, L41, n41(1), n41(2), neigh_41, p1, p2, p3, p4);


% PLOT ELEMENTS %     % PLOT ONE ELEMENT AT A TIME - FOR DEBUGGING 
% PLOT ELEMENTS %     % figure
% PLOT ELEMENTS %     hold on
% PLOT ELEMENTS %   
% PLOT ELEMENTS %   
% PLOT ELEMENTS % % Neighbors $    if neigh_12 > 0
% PLOT ELEMENTS % % Neighbors $      p1n = ele_list_quad(neigh_12,1);
% PLOT ELEMENTS % % Neighbors $      p2n = ele_list_quad(neigh_12,2);
% PLOT ELEMENTS % % Neighbors $      p3n = ele_list_quad(neigh_12,3);
% PLOT ELEMENTS % % Neighbors $      p4n = ele_list_quad(neigh_12,4);
% PLOT ELEMENTS % % Neighbors $  
% PLOT ELEMENTS % % Neighbors $      plot([xP(p1n),xP(p2n),xP(p3n),xP(p4n),xP(p1n)],[yP(p1n),yP(p2n),yP(p3n),yP(p4n),yP(p1n)],'-r','linewidth',2)
% PLOT ELEMENTS % % Neighbors $    end
% PLOT ELEMENTS % % Neighbors $  
% PLOT ELEMENTS % % Neighbors $    if neigh_23 > 0
% PLOT ELEMENTS % % Neighbors $      p1n = ele_list_quad(neigh_23,1);
% PLOT ELEMENTS % % Neighbors $      p2n = ele_list_quad(neigh_23,2);
% PLOT ELEMENTS % % Neighbors $      p3n = ele_list_quad(neigh_23,3);
% PLOT ELEMENTS % % Neighbors $      p4n = ele_list_quad(neigh_23,4);
% PLOT ELEMENTS % % Neighbors $  
% PLOT ELEMENTS % % Neighbors $      plot([xP(p1n),xP(p2n),xP(p3n),xP(p4n),xP(p1n)],[yP(p1n),yP(p2n),yP(p3n),yP(p4n),yP(p1n)],'-r','linewidth',2)
% PLOT ELEMENTS % % Neighbors $    end
% PLOT ELEMENTS % % Neighbors $  
% PLOT ELEMENTS % % Neighbors $    if neigh_34 > 0
% PLOT ELEMENTS % % Neighbors $      p1n = ele_list_quad(neigh_34,1);
% PLOT ELEMENTS % % Neighbors $      p2n = ele_list_quad(neigh_34,2);
% PLOT ELEMENTS % % Neighbors $      p3n = ele_list_quad(neigh_34,3);
% PLOT ELEMENTS % % Neighbors $      p4n = ele_list_quad(neigh_34,4);
% PLOT ELEMENTS % % Neighbors $  
% PLOT ELEMENTS % % Neighbors $      plot([xP(p1n),xP(p2n),xP(p3n),xP(p4n),xP(p1n)],[yP(p1n),yP(p2n),yP(p3n),yP(p4n),yP(p1n)],'-r','linewidth',2)
% PLOT ELEMENTS % % Neighbors $    end
% PLOT ELEMENTS % % Neighbors $  
% PLOT ELEMENTS % % Neighbors $    if neigh_41 > 0
% PLOT ELEMENTS % % Neighbors $      p1n = ele_list_quad(neigh_41,1);
% PLOT ELEMENTS % % Neighbors $      p2n = ele_list_quad(neigh_41,2);
% PLOT ELEMENTS % % Neighbors $      p3n = ele_list_quad(neigh_41,3);
% PLOT ELEMENTS % % Neighbors $      p4n = ele_list_quad(neigh_41,4);
% PLOT ELEMENTS % % Neighbors $  
% PLOT ELEMENTS % % Neighbors $      plot([xP(p1n),xP(p2n),xP(p3n),xP(p4n),xP(p1n)],[yP(p1n),yP(p2n),yP(p3n),yP(p4n),yP(p1n)],'-r','linewidth',2)
% PLOT ELEMENTS % % Neighbors $    end
% PLOT ELEMENTS %  
% PLOT ELEMENTS %     plot([xP(p1),xP(p2),xP(p3),xP(p4),xP(p1)],[yP(p1),yP(p2),yP(p3),yP(p4),yP(p1)],'-k','linewidth',2)
% PLOT ELEMENTS %   
% PLOT ELEMENTS %     if (neigh_12 == BC_inlet_ID)
% PLOT ELEMENTS %       plot([xP(p1),xP(p2)],[yP(p1),yP(p2)],'--b','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_12 == BC_outlet_ID)
% PLOT ELEMENTS %       plot([xP(p1),xP(p2)],[yP(p1),yP(p2)],'--c','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_12 == BC_wall_ID)
% PLOT ELEMENTS %       plot([xP(p1),xP(p2)],[yP(p1),yP(p2)],'--m','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_12 == BC_sym_ID)
% PLOT ELEMENTS %       plot([xP(p1),xP(p2)],[yP(p1),yP(p2)],'--y','linewidth',2)
% PLOT ELEMENTS %   
% PLOT ELEMENTS %       % pause()
% PLOT ELEMENTS %     end
% PLOT ELEMENTS %   
% PLOT ELEMENTS %     if (neigh_23 == BC_inlet_ID)
% PLOT ELEMENTS %       plot([xP(p2),xP(p3)],[yP(p2),yP(p3)],'--b','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_23 == BC_outlet_ID)
% PLOT ELEMENTS %       plot([xP(p2),xP(p3)],[yP(p2),yP(p3)],'--c','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_23 == BC_wall_ID)
% PLOT ELEMENTS %       plot([xP(p2),xP(p3)],[yP(p2),yP(p3)],'--m','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_23 == BC_sym_ID)
% PLOT ELEMENTS %       plot([xP(p2),xP(p3)],[yP(p2),yP(p3)],'--y','linewidth',2)
% PLOT ELEMENTS %   
% PLOT ELEMENTS %       % pause()
% PLOT ELEMENTS %     end
% PLOT ELEMENTS %   
% PLOT ELEMENTS %     if (neigh_34 == BC_inlet_ID)
% PLOT ELEMENTS %       plot([xP(p3),xP(p4)],[yP(p3),yP(p4)],'--b','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_34 == BC_outlet_ID)
% PLOT ELEMENTS %       plot([xP(p3),xP(p4)],[yP(p3),yP(p4)],'--c','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_34 == BC_wall_ID)
% PLOT ELEMENTS %       plot([xP(p3),xP(p4)],[yP(p3),yP(p4)],'--m','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_34 == BC_sym_ID)
% PLOT ELEMENTS %       plot([xP(p3),xP(p4)],[yP(p3),yP(p4)],'--y','linewidth',2)
% PLOT ELEMENTS %   
% PLOT ELEMENTS %       % pause()
% PLOT ELEMENTS %     end
% PLOT ELEMENTS %   
% PLOT ELEMENTS %     if (neigh_41 == BC_inlet_ID)
% PLOT ELEMENTS %       plot([xP(p4),xP(p1)],[yP(p4),yP(p1)],'--b','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_41 == BC_outlet_ID)
% PLOT ELEMENTS %       plot([xP(p4),xP(p1)],[yP(p4),yP(p1)],'--c','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_41 == BC_wall_ID)
% PLOT ELEMENTS %       plot([xP(p4),xP(p1)],[yP(p4),yP(p1)],'--m','linewidth',2)
% PLOT ELEMENTS %     elseif (neigh_41 == BC_sym_ID)
% PLOT ELEMENTS %       plot([xP(p4),xP(p1)],[yP(p4),yP(p1)],'--y','linewidth',2)
% PLOT ELEMENTS %   
% PLOT ELEMENTS %       % pause()
% PLOT ELEMENTS %     end
% PLOT ELEMENTS % %  
% PLOT ELEMENTS % %  
% PLOT ELEMENTS % %    plot(xC,yC,'ok','linewidth',2)
% PLOT ELEMENTS % %  
% PLOT ELEMENTS % %    L12
% PLOT ELEMENTS % %    L23
% PLOT ELEMENTS % %    L34
% PLOT ELEMENTS % %    L41
% PLOT ELEMENTS % %    A
% PLOT ELEMENTS % %  
% PLOT ELEMENTS % %    x12 = (xP(p2) + xP(p1))/2;
% PLOT ELEMENTS % %    y12 = (yP(p2) + yP(p1))/2;
% PLOT ELEMENTS % %    quiver(x12,y12,n12(1)*0.01,n12(2)*0.01)
% PLOT ELEMENTS % %   
% PLOT ELEMENTS % %    x23 = (xP(p3) + xP(p2))/2;
% PLOT ELEMENTS % %    y23 = (yP(p3) + yP(p2))/2;
% PLOT ELEMENTS % %    quiver(x23,y23,n23(1)*0.01,n23(2)*0.01) 
% PLOT ELEMENTS % %  
% PLOT ELEMENTS % %    x34 = (xP(p4) + xP(p3))/2;
% PLOT ELEMENTS % %    y34 = (yP(p4) + yP(p3))/2;
% PLOT ELEMENTS % %    quiver(x34,y34,n34(1)*0.01,n34(2)*0.01) 
% PLOT ELEMENTS % % 
% PLOT ELEMENTS % %    x41 = (xP(p4) + xP(p1))/2;
% PLOT ELEMENTS % %    y41 = (yP(p4) + yP(p1))/2;
% PLOT ELEMENTS % %    quiver(x41,y41,n41(1)*0.01,n41(2)*0.01) 
% PLOT ELEMENTS % % 
% PLOT ELEMENTS %     axis square
% PLOT ELEMENTS %     hold on
% PLOT ELEMENTS % %    pause()

end
toc

fclose(fid_write);

% ++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++




% WRITE VTK QUAD DATA 

fid_VTK = fopen('output.vtk', 'w');

fprintf(fid_VTK,'# vtk DataFile Version 2.0\n');
fprintf(fid_VTK,'Quad data\n');
fprintf(fid_VTK,'ASCII\n');
fprintf(fid_VTK,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid_VTK,'\n');

fprintf(fid_VTK,'POINTS %d float\n', Nnodes);
for ii = 1:Nnodes
  fprintf(fid_VTK,'%f %f 0.0\n', xP(ii), yP(ii));
end

fprintf(fid_VTK,'\n');
fprintf(fid_VTK, 'CELLS %d %d\n', Nelem, Nelem*5); % For triangles, use Nelem*4
for ii = 1:Nelem
  p1 = ele_list_quad(ii,1) - 1; % VTK indices start from 0
  p2 = ele_list_quad(ii,2) - 1; % VTK indices start from 0
  p3 = ele_list_quad(ii,3) - 1; % VTK indices start from 0
  p4 = ele_list_quad(ii,4) - 1; % VTK indices start from 0
  fprintf(fid_VTK, '%d %d %d %d %d\n', 4, p1, p2, p3, p4);
end

fprintf(fid_VTK,'\n');
fprintf(fid_VTK,'CELL_TYPES %d\n', Nelem);
for ii = 1:Nelem
  fprintf(fid_VTK,'9\n'); % 5 for triangles. 9 for quad
end

fprintf(fid_VTK,'\n');

fprintf(fid_VTK,'CELL_DATA %d\n', Nelem);

% CREATE A FIELD
fprintf(fid_VTK,'SCALARS %s float 1\n', 'field_1');
fprintf(fid_VTK,'LOOKUP_TABLE default\n');

for ii = 1:Nelem
  p1 = ele_list_quad(ii,1);
  x = xP(p1);
  y = yP(p1);
  fprintf(fid_VTK,'%f\n', sin(pi*x)+sin(2*pi*y)); % 5 for triangles. 9 for quad
end


fclose(fid_VTK)
% 
% CELL_DATA 2
% SCALARS field1 float 1
% LOOKUP_TABLE default
% 3
% 1
% 
% SCALARS field2 float 1
% LOOKUP_TABLE default
% 1
% 9
% 
