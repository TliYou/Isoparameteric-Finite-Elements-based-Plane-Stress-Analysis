function [globalDispl, globalForce] = solver(globalKmatrix,constrainedDispl, constrainedDispl_globalidx,constrainedForce, totalDOF)

% get reduced matrix      
reduced_Force = constrainedForce;
reduced_Force(constrainedDispl_globalidx,:) = [];

reduced_K = globalKmatrix;
reduced_K(constrainedDispl_globalidx,:) = [];
nozero = find(constrainedDispl);
reduced_Force = reduced_Force - reduced_K(:,nozero)*constrainedDispl(nozero);
reduced_K(:,constrainedDispl_globalidx) = [];
freeDispl = reduced_K \ reduced_Force;

% get global displacement and force
globalDispl = constrainedDispl;
freeDispl_globalidx = setdiff((1:totalDOF), constrainedDispl_globalidx);
globalDispl(freeDispl_globalidx) = freeDispl;
globalForce = globalKmatrix*globalDispl;




