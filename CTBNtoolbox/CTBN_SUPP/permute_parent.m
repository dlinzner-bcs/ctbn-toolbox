function [node]=permute_parent(node0,i,p)

L=length(node0);
node=node0;
node(i).parents=[];

m=randi(L,1,p);
if p>1
    while ((sum(m==i)==1)||sum(unique(m)~=m)~=p)
        m=randi(L,1,p);
    end
else
    while (sum(m==i)==1)
        m=randi(L,1,p);
    end
    
end
node(i).parents=m;


for i=1:L
    node(i).children=[];
end

for i=1:L
    for j=1:length(node(i).parents)
        node(node(i).parents(j)).children=[node(node(i).parents(j)).children, i];
    end
end
end



