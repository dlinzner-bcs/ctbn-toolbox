function [i] = statePos(n,D,L)
%return state position in binary basis
%   Detailed explanation goes here

s=de2bi(0:(D^L)-1);

i=0;
for j=1:D^L
    
    if(s(j,:)-n)==0
    i=j;
    end

end

