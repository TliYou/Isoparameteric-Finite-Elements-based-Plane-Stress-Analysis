clear all; clc; close all

L = 1;       %Length of bar
w = 0.12;    %depth
h = 0.1; 
NU = 0.3;    %poison ratio

P = -3000;   % Applied force on node 22
E = 200e9;

numEle = 2;
numNode = 13;
dof = numNode*2;

%% Define each node coordinate
nodeCoordinate = zeros(numNode,2);

for k=1:numNode
    if k<=5
        nodeCoordinate(k,:) = [0+(k-1)*L/4, 0];
    elseif k>5 && k<=8
        nodeCoordinate(k,:) = [0+(k-6)*L/2, h/2];
    else
        nodeCoordinate(k,:) = [0+(k-9)*L/4, h];
    end
end

%% Define Connectivity

connectivity = [1 3 11 9 2 7 10 6;
                3 5 13 11 4 8 12 7];

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
    x5 = nodeCoordinate(nodes(5),1);
    y5 = nodeCoordinate(nodes(5),2);
    x6 = nodeCoordinate(nodes(6),1);
    y6 = nodeCoordinate(nodes(6),2);
    x7 = nodeCoordinate(nodes(7),1);
    y7 = nodeCoordinate(nodes(7),2);
    x8 = nodeCoordinate(nodes(8),1);
    y8 = nodeCoordinate(nodes(8),2);
    
    k = Q8ElementStiffness(E,NU,w,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8);
    
    K = K + Q8Assemble(K_zero, k, nodes);
end

%% Boundary and Compatibility Conditions 

F = zeros(dof,1);
F(13*2) = P;

d = zeros(dof,1);

c1 = 1; c2 = 6; c3 =9;

nxt = @(c)(2*c+1);     %Next nodal location 
bfr = @(c) (2*(c-1));  %Before nodal location 

 
% K_mod = [K((c1+1)*2-1:(c2-1)*2,(c1+1)*2-1:(c2-1)*2) ,   K((c1+1)*2-1:(c2-1)*2, (c2+1)*2-1:(c3-1)*2), K((c1+1)*2-1:(c2-1)*2, (c3+1)*2-1:dof);
%          K((c2+1)*2-1:dof, (c1+1)*2-1:(c2-1)*2),             K((c2+1)*2-1:dof, (c2+1)*2-1:dof)];

K_mod =  [K(nxt(c1):bfr(c2), nxt(c1):bfr(c2)), K(nxt(c1):bfr(c2), nxt(c2):bfr(c3)), K(nxt(c1):bfr(c2), nxt(c3):dof);
          K(nxt(c2):bfr(c3), nxt(c1):bfr(c2)), K(nxt(c2):bfr(c3), nxt(c2):bfr(c3)), K(nxt(c2):bfr(c3), nxt(c3):dof);
          K(nxt(c3):dof,     nxt(c1):bfr(c2)), K(nxt(c3):dof,     nxt(c2):bfr(c3)), K(nxt(c3):dof,     nxt(c3):dof)];


F_mod = [F(nxt(c1):bfr(c2)); F(nxt(c2):bfr(c3)); F(nxt(c3):dof)];

d_mod = K_mod\F_mod;
% 
fprintf("Vertical displacement at node 13: %.10f m\n\n", d_mod(end))

d = [d(1:2);
     d_mod(1:8)
     d(11:12);
     d_mod(9:12)
     d(18:19);
     d_mod(13:end)];
 
 F = K*d;
 
fprintf('node\t x-dis[m]\t y-dis[m]\n')
 for no=1:numNode
    fprintf('%d   	 %.4e    %.4e\n',no, d(2*no-1), d(2*no))
 end
 
  fprintf('\n\nnode\t F_x[N]\t\t  F_y[N]\n')
 for no=1:numNode
    fprintf('%d   	 %.4e    %.4e\n',no, F(2*no-1), F(2*no))
 end





    
    
    






