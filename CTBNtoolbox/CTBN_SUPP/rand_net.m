function [B] = rand_net(L,MAX_PAR)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

B=zeros(L);

for l=1:L
    m=[1:L];
    par=randi(MAX_PAR+1)-1;   
    for i=1:par    
        j=randi(L-(i-1));
        B(l,m(j))=1;
        m(j)=[]; 
    end 
end

end

