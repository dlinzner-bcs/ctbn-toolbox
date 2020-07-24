function sample=sampleOfOneSpin(spin,steps,numbOfV,node, timeInterval, stepSize)
%input: spin - number of spin (node)
%       steps - amount of steps in simulation     
%       numbOfV - number of vertices (number of nodes)
%       node - lib of nodes
%       timeInterval - integer that means end of time interval of
%       simulation (the bigger number of nodes is the smaller it shoul be
%       compareing to amount of steps)
%       stepSize - accuracy
%output: vector of spin`s (node`s) sates

[timeVector, statesMatrix]=changesInNetwork(steps,numbOfV,node);
matrixForMagnetization=createMatrixForM(timeVector, statesMatrix, timeInterval, stepSize);
for w=1:timeInterval/stepSize
    sample(w)=matrixForMagnetization(w,spin);
end
end