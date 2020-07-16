clear; close all; clc;

%% Set input data
[inputData]=EulerBernouli;

%% Get element stiffness data
[elementData] = getElementStiffness(inputData);

%% Assemble global stiffness matrix
[globalKmatrix] = assemGlobalStiffness(elementData, inputData);

%% Reduce matrix and solve the equation
[outputData] = solver(globalKmatrix, inputData);

%% Postprocessing
postProcessing(inputData, elementData, outputData);

fprintf('\n');

