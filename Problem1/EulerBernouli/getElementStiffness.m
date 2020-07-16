function [elementData] = getElementStiffness(inputData)

numberofElement = inputData.numberofElement;
nodeConnectivity = inputData.nodeConnectivity;
DOFperNode = inputData.DOFperNode;
nodeCoordinate = inputData.nodeCoordinate;
XsectionArea = inputData.XsectionArea;
elasticModulus = inputData.elasticModulus;
momentofInertia = inputData.momentofInertia;

for en = 1:numberofElement
    % Generate global stiffness matrix of element
    elementNodes = nodeConnectivity(en,:);
    positionVector1 = nodeCoordinate(elementNodes(1),:);
    positionVector2 = nodeCoordinate(elementNodes(2),:);
    L = norm(positionVector2 - positionVector1);
    unitVector = (positionVector2 - positionVector1)/L;
    
    Area = XsectionArea(en);
    constitutiveMat = elasticModulus(en);
    C1 = constitutiveMat*Area/L;
    C2 = (momentofInertia(en)*constitutiveMat)./L^3;
    
    localKmatrix = [C1   0       0       -C1      0     0;
                     0   12*C2   6*L*C2   0  -12*C2   6*L*C2;
                     0    6*L*C2 4*L^2*C2 0  -6*L*C2  2*L^2*C2;
                     -C1  0      0        C1  0       0;
                     0   -12*C2 -6*L*C2    0  12*C2  -6*L*C2;
                     0    6*L*C2 2*L^2*C2  0 -6*L*C2 4*L^2*C2];
    
    e_x = [1 0];
    e_y = [0 1];
    C = dot(unitVector, e_x);
    S = dot(unitVector, e_y);
    
    transformationMat = [C S 0 0 0 0;
                        -S C 0 0 0 0;
                         0 0 1 0 0 0;
                         0 0 0 C S 0;
                         0 0 0 -S C 0;
                         0 0 0 0 0 1];
        
    globalKmatrix_element = transformationMat' * localKmatrix * transformationMat;

    % Index matrix with global index
    global_idx = [];
    
    for i = 1:length(elementNodes)
        node = elementNodes(i);
        global_idx = [ global_idx DOFperNode*(node-1)+1 : DOFperNode*node ];
    end
    globalKdata(en,:) = { global_idx, globalKmatrix_element };
    localKdata(en,:) = { elementNodes, localKmatrix, transformationMat,...
        constitutiveMat, L};
end

elementData.globalKdata = globalKdata;
elementData.localKdata = localKdata;

