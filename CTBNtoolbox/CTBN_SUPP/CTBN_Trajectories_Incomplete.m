function [time0,DATA0] = CTBN_Trajectories_Incomplete(node0,steps,steps0,dt0)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
L=length(node0);

[timeVector, statesMatrix]=changesInNetwork(steps,L,node0);
[DATA,time]=createMatrixForM(timeVector, statesMatrix, timeVector(end)-timeVector(1), dt0);

sm=size(DATA);

MM=[1:sm];
MM1=MM(unique(randi(sm(1),steps0,1)));
DATA0=DATA(MM1,:);
time0=time(MM1);

end

