clear all; clc; close all
% Dimensions & Properties
L = 56;    %[mm]
t=1;       %[mm]
E = 7e4;   %[N/mm^2]
NU = 0.25;
D = 20;    %[mm]

% Quarter Model rectangle vertex locations
x = [0 L/2  L/2 0];
y = [0  0   D/2 D/2];

% rectangle
rect = [3 4 x y]';

% ration: rho/d
ratio = 0.1:0.1:1.0;

%% Plot of Meshes onto the geometries

for lp = 1:length(ratio)
    
    d = D/(2*ratio(lp)+1);
    r = (D-d)/2;
    fprintf('Check the ratio: r/d = %.2f\n', r/d);
    
    
    % Circle [1 x_c y_c r]
    C = [1 0 0 r]';
    % Zero Padding
    C = [C; zeros(length(rect)-length(C),1)];
    
    % Matrix with all geometry description
    gd = [rect C];
    
    % Create name space
    ns = char('rect', 'C');
    ns = ns';
    
    % Set Geometry Formula
    sf = 'rect - C';
    
    [dl, bt] = decsg(gd,sf,ns);
    
    % Use PDE plotting tools
    model = createpde;
    geometryFromEdges(model,dl);
    
    mesh_mod = generateMesh(model, 'Hmax', 2.5, 'Hmin', 0.1, 'GeometricOrder','linear');
    meshAll(lp)= mesh_mod;            %Save meshes for all cases for later analysis
%     figure(lp)
%     pdemesh(mesh_mod)
%     axis([-4 L/2+4 -4 D/2+4])
%     title(['Mesh plot for \rho/d = ',num2str(ratio(lp))])
end

%% Identification of boundary nodes
for lp=1:length(ratio)
    d = D/(2*ratio(lp)+1);
    r = (D-d)/2;
    
    nodes_lower = findNodes(meshAll(lp), 'box', [r L/2], [-0.1 0.1]);   %nodes on right side
    nodes_left = findNodes(meshAll(lp), 'box', [0 0.1], [r D/2]);       %nodes on left side
    nodes_right = findNodes(meshAll(lp),'box',[L/2-0.1,L/2+0.1],[0;D/2]);  %nodes on right side

    figure(lp)
    pdemesh(meshAll(lp))
    hold on
    plot(meshAll(lp).Nodes(1,nodes_lower),meshAll(lp).Nodes(2,nodes_lower),'or','MarkerFaceColor','g')
    for go=1:length(nodes_lower)
        text(meshAll(lp).Nodes(1,nodes_lower(go))+0.4,meshAll(lp).Nodes(2,nodes_lower(go)), num2str(nodes_lower(go)));
    end
    
    plot(meshAll(lp).Nodes(1,nodes_left), meshAll(lp).Nodes(2,nodes_left),'or','MarkerFaceColor','b')
    for go=1:length(nodes_left)
        text(meshAll(lp).Nodes(1,nodes_left(go))+0.4,meshAll(lp).Nodes(2,nodes_left(go)), num2str(nodes_left(go)));
    end
    
    plot(meshAll(lp).Nodes(1,nodes_right),   meshAll(lp).Nodes(2,nodes_right),'or','MarkerFaceColor','m')
    for go=1:length(nodes_right)
        text(meshAll(lp).Nodes(1,nodes_right(go))+0.4,meshAll(lp).Nodes(2,nodes_right(go)), num2str(nodes_right(go)));
    end
    
    axis([-4 L/2+4 -4 D/2+4])
    title(['Boundary nodes for \rho/d = ',num2str(ratio(lp))])
    axis equal
    
    %Element and Global Stiffness Matrix
    numEle = size(meshAll(lp).Elements,2);
    numNode = size(meshAll(lp).Nodes,2);
    dof = numNode*2;
    
    nodeCoordinate = meshAll(lp).Nodes';
    connectivity = meshAll(lp).Elements';
    
    K_zero = zeros(dof,dof);
    K = zeros(dof,dof);
    
    for g = 1:numEle
        nodes = connectivity(g,:);
        xi = nodeCoordinate(nodes(1),1);
        yi = nodeCoordinate(nodes(1),2);
        xj = nodeCoordinate(nodes(2),1);
        yj = nodeCoordinate(nodes(2),2);
        xm = nodeCoordinate(nodes(3),1);
        ym = nodeCoordinate(nodes(3),2);
        k = CSTElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym);
        
        K = K + CSTAssemble(K_zero, k, nodes);
        k_element(lp,g) = {k};  % saves the element striffness matrices
    end
    K_global(lp) = {K};         %saves the global stiffness matrices for each case
    
end

% Global Displacement Calculation
for lp=1:length(ratio)
    d = D/(2*ratio(lp)+1);
    r = (D-d)/2;
    
    nodes_lower = findNodes(meshAll(lp), 'box', [r L/2], [-0.1 0.1]);      %nodes on right side
    nodes_left = findNodes(meshAll(lp), 'box', [0 0.1], [r D/2]);          %nodes on left side
    nodes_right = findNodes(meshAll(lp),'box',[L/2-0.1,L/2+0.1],[0;D/2]);  %nodes on right side
    
    numEle = size(meshAll(lp).Elements,2);
    numNode = size(meshAll(lp).Nodes,2);
    dof = numNode*2;
    
    % Boundary and Compatibility Conditions
    constrainedForce = zeros(dof,1);
    constrainedDispl = zeros(dof,1);
    constrainedDispl_globalidx = [2*nodes_left-1 2*nodes_lower];
    
    F_e = 31.25 ; F_m = 2*31.25;
    constrainedForce(2*nodes_right'-1) = [F_e F_e F_m*ones(1,length(nodes_right)-2)]';
    
    K = cell2mat(K_global(lp));
    
    [globalDispl, globalForce] = solver(K,constrainedDispl, constrainedDispl_globalidx,constrainedForce, dof);
    globalDispl_all(lp) = {globalDispl};   %Saves the global displacement data for all cases
end

% Calculation of stress in each element
for lp = 1:length(ratio)
    globalDispl = cell2mat(globalDispl_all(lp));
    
    numEle = size(meshAll(lp).Elements,2);
    numNode = size(meshAll(lp).Nodes,2);
    dof = numNode*2;
    
    nodeCoordinate = meshAll(lp).Nodes';
    connectivity = meshAll(lp).Elements';
    
    stress_mode = [];
    for g = 1:numEle
        k = cell2mat(k_element(lp,g));
        nodes = connectivity(g,:);
        d = globalDispl([2*nodes(1)-1, 2*nodes(1), 2*nodes(2)-1, 2*nodes(2), 2*nodes(3)-1, 2*nodes(3)]');  %local dispalcement d for each element
        
        xi = nodeCoordinate(nodes(1),1);
        yi = nodeCoordinate(nodes(1),2);
        xj = nodeCoordinate(nodes(2),1);
        yj = nodeCoordinate(nodes(2),2);
        xm = nodeCoordinate(nodes(3),1);
        ym = nodeCoordinate(nodes(3),2);
        stress = CSTElementStresses(E,NU,xi,yi,xj,yj,xm,ym,d);
        stress_mode = [stress_mode, stress];
    end
    stress_all(lp) = {stress_mode};
end


% Determinination of maximum normal stress and the nominal stress

max_stress = [];
stress_avg = [];
Scc = [];
fprintf('\nr/d \t Max Normal stress[N/mm^2]\t Nominal Stress[N/mm^2]\t\t Scc\n');
fprintf('=======================================================================\n')

for lp=1:length(ratio)
    d = D/(2*ratio(lp)+1);
    stress_a = 500/d*1;
    stress = cell2mat(stress_all(lp));
    max_s = max(abs(stress(1,:)));               %maximum stress in the x-direction
    Scc_i = max_s/stress_a;                      %Stress Concentration factor
    fprintf('%.2f\t\t   %.4e\t\t\t       %.4e\t\t\t  %.3f\n', ratio(lp), max_s, stress_a, Scc_i);
    max_stress = [max_stress, max_s];
    stress_avg = [stress_avg, stress_a];
    Scc = [Scc, Scc_i];
end

figure
plot(ratio, Scc,'--o', 'LineWidth',1)
xlabel('Ratio, \rho/d');
ylabel('Stress Concentration Factor S_{cc}')










