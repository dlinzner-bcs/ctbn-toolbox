function [edges] = node2net(node)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
L=length(node);

edges=zeros(L);

for i=1:L
    
    edges(i,node(i).parents)=1;
     edges(node(i).children,i)=1;
    
end


end

