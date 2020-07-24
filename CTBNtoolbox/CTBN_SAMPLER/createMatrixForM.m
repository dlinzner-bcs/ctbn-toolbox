function [matrixForMagnetization,time]=createMatrixForM(timeVector, statesMatrix, timeInterval, intervalOfTimeGap)
%creates new full matrix of states of the network on time interval that
%needed (should be bigger than amount of steps thak network was changing) 
lim=length(timeVector);
elemOfM=1;
timeT=0;

for w=1:timeInterval/intervalOfTimeGap
     time(w)=timeT;
    if elemOfM>=lim
        matrixForMagnetization(w,:)=statesMatrix(lim,:);
    else
        if timeT<timeVector(elemOfM)
            matrixForMagnetization(w,:)=statesMatrix(elemOfM,:);
            
        else
            matrixForMagnetization(w,:)=statesMatrix(elemOfM+1,:);
            elemOfM=elemOfM+1;
        end
        timeT=timeT+intervalOfTimeGap;
    end
   
end
end