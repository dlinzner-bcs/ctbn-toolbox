function parentsState = checkParentsState(node, i)
%returns a vector of parent`s state of node i
n=length(node(i).parents);
parentsState=0;
if n>0
    parentsState=node(node(i).parents(1)).state;
    for j=2:n
        r=node(node(i).parents(j)).state;
        parentsState=[parentsState,r];
    end
end


end