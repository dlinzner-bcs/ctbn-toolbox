function [EZ,X0] = CTBN2CTMC_data(statesMatrix)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 sz=size(statesMatrix);
 L=sz(3);
    s=de2bi(0:(2^L)-1);
    EZ=zeros(sz(2),2^L);
    for n=1:sz(2)
        for j=1:length(s)
            
            E=zeros(1,2^L);
            E(statePos(s(j,:),2,L))=1;
            a=1;
            
            for i=1:L
                a=a*statesMatrix(s(j,i)+1,n,i);
            end
            EZ(n,:)=EZ(n,:)+a*E;
        end
    end
    
    X0=(1/(2^L))*ones(1,2^L);

end

