function [node]=createLibOfNodesGammaPrior(L,A,init,tau,M)
%input: L - number of vertices (number of nodes)
%       edges - adjency matrix (node i is a child of node j if edges(i,j)=1)
%       init - vector of states which consist of 1 and -1 only relevant for sampling
%       tau - prior dwelling-times
%       M - prior transition number
%output: library of structures node{i} with random rate matrices, drawn
%from gamma prior with hyperparams tau and M

%empty lib of nodes
lib(L,'node')=struct('state',0);

% for each node create families
for i=1:L
    node(i).state=init(i);
    node(i).parents=[];
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
%create conditional rate matrixes
for w=1:L
    if length(node(w).parents)>0
        for v=1:2^length(node(w).parents)
            for u=1:length(node(w).parents)
                if node(w).allPosStatesOfParents(v,u)==2
                    node(w).allPosStatesOfParents(v,u)=-1;
                end           
            end
            wDown = gamrnd(M(v,w,1,2),1/tau(v,w,2));
            wUp   = gamrnd(M(v,w,2,1),1/tau(v,w,1));
            cellOfCRM{v}=[-wDown, wDown; wUp, -wUp];
        end
    else
        wDown = gamrnd(M(1,w,1,2),1/tau(1,w,2));
        wUp   = gamrnd(M(1,w,2,1),1/tau(1,w,1));
        cellOfCRM{1}=[-wDown, wDown; wUp, -wUp];
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