function node=createLibOfNodes(L,A, init, temp,tau)
%input: L - number of vertices (number of nodes)
%       A - adjency matrix (node i is a child of node j if edges(i,j)=1)
%       init - vector of states which consist of 1 and -1
%       temp - temperature
%output: library of structures node{i} with rate matrices encoding glauber
%dynamics

%empty lib of nodes
lib(L,'node')=struct('state',0);

for i=1:L
    node(i).parents=[];
end

% for each node create families
for i=1:L    
    node(i).state=init(i);
    %vector of parents
    connection=A(i,:);
    firstParent=true;
    for j=1:L
        if connection(j)==1 & firstParent==false
            a=node(i).parents;
            node(i).parents=[a, j];
        end
        if firstParent & connection(j)==1
            node(i).parents=[j];
            firstParent=false;
        end
    end    
end

for i=1:L
    %matrix of all possible states of parents
    if isempty(node(i).parents)  
        node(i).allPosStatesOfParents=[];
    else
        parentsLength=length(node(i).parents);
        if parentsLength>=1
            s=de2bi(0:(2^parentsLength)-1);
            s=s+1;
            node(i).allPosStatesOfParents=s;
        else
            node(i).allPosStatesOfParents=0;
        end
    end
end

%creating conditional rate matrixes
for w=1:L
    for v=1:2^length(node(w).parents)
        if length(node(w).parents)>0
            for u=1:length(node(w).parents)
                if node(w).allPosStatesOfParents(v,u)==2
                    node(w).allPosStatesOfParents(v,u)=-1;
                end
                sumOfParSt=sum(node(w).allPosStatesOfParents(v,:));                
            end
        else
            sumOfParSt=0;
        end       
        wDown=(1/2)*(1+tanh(temp*sumOfParSt));
        wUp=(1/2)*(1-tanh(temp*sumOfParSt));
        cellOfCRM{v}=tau*[-wDown, wDown; wUp, -wUp];    
    end   
    node(w).cellOfCondRM=cellOfCRM;
    clear cellOfCRM;
end

%create children consistent with parents
for i=1:L
    node(i).children=[];
end

for i=1:L
    for j=1:length(node(i).parents)
        node(node(i).parents(j)).children=[node(node(i).parents(j)).children, i];
    end
end
end