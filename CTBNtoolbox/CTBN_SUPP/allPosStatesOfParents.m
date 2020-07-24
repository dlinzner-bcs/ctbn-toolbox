function [node] = allPosStatesOfParents(node)
%input: node - NonParamCTBN
%       i    - index of initialized node

%output: All possible states of parents A
for n=1:numel(node)
    
    parentsLength=length(node(n).parents);
    
    %create matrix of all possible states of parents
    U=[1];
    for j=node(n).parents
        
        B=[1:length(node(j).Omega)];
        U=combvec(U,B);
        
    end
    
    
    if parentsLength>=1
        U(1,:)=[];
    else
        U=[];
    end
    
    node(n).allPosStatesOfParents=U';
    
end