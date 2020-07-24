function choosedRateM= chooseConditionalRM(parentsState, node, i)
%returns conditional rate matrix depending of node`s parent`s state
choosedRateM=0;

for j=1:length(node(i).allPosStatesOfParents)
    if parentsState==node(i).allPosStatesOfParents(j,:)
        choosedRateM=node(i).cellOfCondRM{j};
    end

end