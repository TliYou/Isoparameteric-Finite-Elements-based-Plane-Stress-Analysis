function postProcessing(inputData, elementData, outputData)
 
totalDOF = inputData.totalDOF;
numberofElement = inputData.numberofElement;
globalDispl = outputData.globalDispl;
nodeConnectivity = inputData.nodeConnectivity;
nodeCoordinate = inputData.nodeCoordinate;
localKdata = elementData.localKdata;
globalKdata= elementData.globalKdata;

globalForce = outputData.globalForce;

% Display configuration

for en = 1:numberofElement
    elementNodes = nodeConnectivity(en,:);
    xx = nodeCoordinate(elementNodes,1);
    yy = nodeCoordinate(elementNodes,2);
    plot(xx,yy,'--','Color','black','MarkerFaceColor','r','MarkerSize',1,'LineWidth',1);
    hold on
end


% Rearrange the global displacement matrix for better manupulation
for i = 1:length(globalDispl)/3
    globalDispl_reshape(i,1) = globalDispl(3*i-2);
    globalDispl_reshape(i,2) = globalDispl(3*i-1);
    globalDispl_reshape(i,3) = globalDispl(3*i);
end

scaleFactor = 20;
for en = 1:numberofElement
    elementNodes = nodeConnectivity(en,:);
    xx = nodeCoordinate(elementNodes,1) + globalDispl_reshape(elementNodes,1)*scaleFactor;
    yy = nodeCoordinate(elementNodes,2) + globalDispl_reshape(elementNodes,2)*scaleFactor;
    plot(xx,yy,'-o','Color','black','MarkerFaceColor','r','MarkerSize',5,'LineWidth',1);
    hold on
end
axis equal
grid on
xlabel('x axis [m]')
ylabel('y axis [m]')
title('Deformation Scale: x20')

fprintf("Vertical displacement at node 2: %.10f m\n\n", globalDispl(5))

fprintf('node\t x-dis[m]\t y-dis[m]\t Z-rot[rad]\n')
 
 for no=1:2
    fprintf('%d   	 %.4e  %.4e  %.4e\n',no, globalDispl(3*no-2), globalDispl(3*no-1), globalDispl(3*no));
 end
 
  fprintf('\n\nnode\t F_x[N]\t\t F_y[N]\t\t\tT_z[Nm]\n')
 for no=1:2
    fprintf('%d   	 %.4e	 %.4e \t%.4e\n',no, globalForce(3*no-2), globalForce(3*no-1),globalForce(3*no))
 end






