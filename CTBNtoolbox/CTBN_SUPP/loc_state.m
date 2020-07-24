function [ind] = loc_state(i,m,D,L)
%return state position in binary basis
%   Detailed explanation goes here

s=de2bi(0:(D^L)-1);

ind=[];

for k=1:D^L
    
    if s(k,i)==m
       
       ind= [ind k];
        
    end
    
end

end
