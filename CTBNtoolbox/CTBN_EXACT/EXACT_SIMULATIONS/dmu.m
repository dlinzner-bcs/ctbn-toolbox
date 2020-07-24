function [ dmu] = dmu(D,L,t,y,ft,f)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

for i=1:D^L
    for j=1:D^L
        
        h(i,j) = interp1(ft, f(:,i,j), t);
    
    end
end

dmu=h*y;

end

