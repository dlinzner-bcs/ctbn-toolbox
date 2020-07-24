function [timeVector, statesMatrix]=changesInNetworkD(steps,numbOfV,node)
%function that makes one simulation of network 
%input: steps - amount of steps in simulation     
%       numbOfV - number of vertices (number of nodes)
%       node - lib of nodes
%output: timeVector contains time when there was changes in network
t=0;
for q=1:steps
    
    [rates, generalRate]=getRatesD(numbOfV,node);
    for s=1:numbOfV
        networkState(s)=node(s).state;
    end
    
    t=t+exprnd(1/generalRate);
    statesMatrix(q,:)=networkState;
    timeVector(q)=t;

    j=chooseNodeToSwitch(rates, generalRate);
    parentsState = checkParentsState(node, j);
    choosedRateM= chooseConditionalRM(parentsState, node, j);
    
    node(j).currentCondRM=choosedRateM;
    rq=[];
    for d=1:length(node(j).Omega)
       if d~=node(j).state
        rq=[rq sum(rq)+node(j).currentCondRM(node(j).state,d)];
       else
        rq=[rq 10^(-18)];   
       end
    
    end
    rq=rq/max(rq);
    d_=find(rand<rq,1);
    node(j).state=d_;
end
end