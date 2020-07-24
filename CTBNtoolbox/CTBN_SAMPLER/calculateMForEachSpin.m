
function magnetization=calculateMForEachSpin(spin,steps,numbOfV,node, timeInterval, numbOfIter, beginningState, stepSize)
%input: spin - number of spin (node)
%       steps - amount of steps in simulation     
%       numbOfV - number of vertices (number of nodes)
%       node - lib of nodes
%       timeInterval - integer that means end of time interval of
%       simulation (the bigger number of nodes is the smaller it shoul be
%       compareing to amount of steps)
%       numbOfIter - nubmer of iterations (amount of simulations)
%       beginning state - vector of states which consist of 1 and -1
%       stepSize - accuracy
%output: vector of magnetization with given step size for given spin (node)
for v=1:numbOfIter
    [timeVector, statesMatrix]=changesInNetwork(steps,numbOfV,node);
    matrixForMagnetization=createMatrixForM(timeVector, statesMatrix, timeInterval, stepSize);
    for w=1:timeInterval/stepSize
        magnetiz(w,v)=matrixForMagnetization(w,spin);
    end
    for i=1:numbOfV
        node(i).state=beginningState(i);
    end
end
for w=1:timeInterval/stepSize
    magnetization(w)=mean(magnetiz(w,:));
end
end
    