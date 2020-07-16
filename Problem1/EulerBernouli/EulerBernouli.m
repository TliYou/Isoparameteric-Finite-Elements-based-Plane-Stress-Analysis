function [inputData]=EulerBernouli
% Set model data
numberofNode = 2;                       % number of node
numberofElement = 1;                    % number of element
DOFperNode = 3;                         % DOF per node
totalDOF = DOFperNode*numberofNode;     % total DOF
nodeConnectivity(1,:) = [1 2];          % connectivity of element 1

nodeCoordinate(1,:) = [0 0.1];          % coordinate of node 1 [m]
nodeCoordinate(2,:) = [1 0.1];          % coordinate of node 2 [m]

b = 0.12; h=0.1;
XsectionArea = b*h*0.12*ones(1,numberofElement);                 % cross-section area of all elements
elasticModulus = 200e9*ones(1,numberofElement);               % elastic modulus of all elements [Pa]
momentofInertia = (b*h^3/12)*ones(1,numberofElement);               % moment inertia of all elements in SI unit


% Set boundary and loading condition
constrainedDispl = zeros(totalDOF,1);                  % constrained displacement [m]
constrainedDispl_globalidx = [1 2 3];         % global index having constrained displacement
constrainedForce = zeros(totalDOF,1);                  % constrained force [N]
constrainedForce(5) = -3000;


% Save input data
inputData.numberofNode = numberofNode;
inputData.numberofElement = numberofElement;
inputData.DOFperNode = DOFperNode;
inputData.totalDOF = totalDOF;
inputData.nodeConnectivity = nodeConnectivity;
inputData.constrainedDispl = constrainedDispl;
inputData.constrainedDispl_globalidx = constrainedDispl_globalidx;
inputData.constrainedForce = constrainedForce;
inputData.nodeCoordinate = nodeCoordinate;
inputData.XsectionArea = XsectionArea;
inputData.elasticModulus = elasticModulus;
inputData.momentofInertia = momentofInertia;
