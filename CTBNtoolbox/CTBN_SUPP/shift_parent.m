function [node]=shift_parent(node0,i,j)
%shift parents of i to j 
%j can be vector
 
L=length(node0);
node=node0;
node(i).parents=[];
node(i).parents=j;


for i=1:L
    node(i).children=[];
end

for i=1:L
    for j=1:length(node(i).parents)
        node(node(i).parents(j)).children=[node(node(i).parents(j)).children, i];
    end
end
end



