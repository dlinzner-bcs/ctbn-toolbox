function [timeVector, statesMatrix]=changesInNetwork(steps,numbOfV,node)
%function that makes one simulation of network 
%input: steps - amount of steps in simulation     
%       numbOfV - number of vertices (number of nodes)
%       node - lib of nodes
%output: timeVector contains time when there was changes in network
t=0;
for q=1:steps
    
    [rates, generalRate]=getRates(numbOfV,node);
    for s=1:numbOfV
        networkState(s)=node(s).state;
    end
    
    t=t+exprnd(1/generalRate);
    statesMatrix(q,:)=networkState;
    timeVector(q)=t;
    
    choosedNodeNumb=chooseNodeToSwitch(rates, generalRate);
    if node(choosedNodeNumb).state==1
        node(choosedNodeNumb).state=-1;
    else 
        node(choosedNodeNumb).state=1;
    end
    
end
end