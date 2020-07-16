clear all; clc; close all

L = 1;       %Length of bar
w = 0.12;    %depth
h = 0.1; 
NU = 0.3;    %poison ratio

P = -3000;   % Applied force on node 22
E = 200e9;

numEle = 10;
numNode = 22;
dof = numNode*2;

%% Define each node coordinate
nodeCoordinate = zeros(numNode,2);

for k=1:numNode
    if k<=11
        nodeCoordinate(k,:) = [0+(k-1)*L/numEle, 0];
    else
        nodeCoordinate(k,:) = [0+(k-12)*L/numEle, h];
    end
end

%% Define Connectivity
lower = 1:11;
upper = 12:22;
connectivity = zeros(numEle,4); 
for k=1:numEle
    connectivity(k,:) = [lower(k:k+1) upper(k+1:-1:k)];
end

%% Find the global stiffness matrix
clear k

K_zero = zeros(dof,dof);
K= zeros(dof,dof);

for g = 1:numEle
    nodes = connectivity(g,:);
    x1 = nodeCoordinate(nodes(1),1);
    y1 = nodeCoordinate(nodes(1),2);
    x2 = nodeCoordinate(nodes(2),1);
    y2 = nodeCoordinate(nodes(2),2);
    x3 = nodeCoordinate(nodes(3),1);
    y3 = nodeCoordinate(nodes(3),2);
    x4 = nodeCoordinate(nodes(4),1);
    y4 = nodeCoordinate(nodes(4),2);
   
    k = Q6ElementStiffness(E,NU,w,x1,y1,x2,y2,x3,y3,x4,y4);
    
    K = K + Q6Assemble(K_zero, k, nodes);
end

%% Boundary and Compatibility Conditions 

F = zeros(dof,1);
F(22*2) = P;

d = zeros(dof,1);

c1 = 1; c2 = 12;
% constrained_globalidx = [2*c1-1 2*c1 2*c2-1 2*c2];

K_mod = [K((c1+1)*2-1:(c2-1)*2,(c1+1)*2-1:(c2-1)*2) ,   K((c1+1)*2-1:(c2-1)*2, (c2+1)*2-1:dof);
         K((c2+1)*2-1:dof, (c1+1)*2-1:(c2-1)*2),             K((c2+1)*2-1:dof, (c2+1)*2-1:dof)];
F_mod = [F((c1+1)*2-1:(c2-1)*2); F((c2+1)*2-1:dof)];

d_mod = K_mod\F_mod;

fprintf("Vertical displacement at node 22: %.10f m\n\n", d_mod(end))

d = [d(c1*2-1:c1*2);
     d_mod(1:20);
     d(c2*2-1:c2*2);
     d_mod(21:40)];
 
 F = K*d;
 
fprintf('node\t x-dis[m]\t y-dis[m]\n')
 for no=1:numNode
    fprintf('%d   	 %.4e    %.4e\n',no, d(2*no-1), d(2*no))
 end
 
  fprintf('\n\nnode\t F_x[N]\t\t  F_y[N]\n')
 for no=1:numNode
    fprintf('%d   	 %.4e    %.4e\n',no, F(2*no-1), F(2*no))
 end






    
    
    






